import Pkg
Pkg.activate("Modron")

using ProgressMeter
# using Combinatorics
using Kmers
using ArgParse

include("Kct.jl")

abstract type Abstract_Node end

mutable struct Tree_Node <: Abstract_Node
    kct_index::UInt32
    to::Vector{Tree_Node}
    samples::Vector{UInt32}
end
get_seq(node::Abstract_Node, kct::Kct{N}) where {N} = String(kct[node.kct_index][1])

mutable struct Tree_Root <: Abstract_Node
    seed::String
    to::Vector{Abstract_Node}
    samples::Vector{Int}
end

get_seq(node::Tree_Root, kct::Kct{N}) where {N} = node.seed

const nucleotides = ('A', 'T', 'G', 'C')

function parse_commands()
    arg_settings = ArgParseSettings()
    @add_arg_table arg_settings begin
        "walk"
            help = "Signals k-mer walking mode."
            action = :command
        "prepare"
            help = "Signals kct-preparation mode."
            action = :command
    end

    @add_arg_table arg_settings["walk"] begin
        "--seed", "-s"
            help = "The seed to extend by k-mer walking."
            arg_type = String
            required = true
        "--library", "-l"
            help = "The KCT library to use for k-mer walking."
            arg_type = String
            required = true
        "--output", "-o"
            help = "Output file. Will be in fasta format."
            arg_type = String
            required = true
        "--max_recursion", "-r"
            help = "How much the seed can be extended by k-mer walking."
            arg_type = Int
            default = 500
        "--min_count", "-c"
            help = "The minimum count for a k-mer to be admissible. Should be >0."
            arg_type = Int
            default = 1
    end

    @add_arg_table arg_settings["prepare"] begin
        "--data", "-d"
            help = "Folder that contains ONLY jellyfish files to turn into a Kct."
            arg_type = String
            required = true
        "--output", "-o"
            help = "The output kct file. Will be in binary format."
            arg_type = String
            required = true
        "--jellyfish", "-j"
            help = "Path of jellyfish executable"
            arg_type = String
            required = true   
        "-k"
            help = "The size of k-mers in jellyfish files"
            arg_type = Int
            default=31
    end   

    return parse_args(arg_settings)
end


function build_kct(jf_files::Vector{String}, output::String, jellyfish_path::String; k::Int64=31)
    if length(jf_files) == 1
        save(backup_kct(jf_files[1], jellyfish_path, k=k), output)
        return
    end

    # Init stuff
    prog = Progress(length(jf_files), "Parsing Jellyfish Files (binary tree version)"); update!(prog)

    # Stacks previously loaded Kct
    tree_stack = Array{Kct, 1}()

    # Using a stack to emulate a binary tree merge. Adding a new leaf at every loop
    for file in jf_files
        push!(tree_stack, backup_kct(file, k=k, jellyfish_path))
        
        # Then checking if that leaf can be merged with another leaf at the same level.
        # Keep merging up until no matching leaf at that level.
        while length(tree_stack) >= 2 && typeof(tree_stack[end]) == typeof(tree_stack[end-1])
            update!(prog, showvalues = [("Tree size on RAM: ", Base.format_bytes(Base.summarysize(tree_stack))),
                                        ("Merging", "["*join([typeof(k) for k in tree_stack[1:end-1]], " | ")*" ← $(typeof(tree_stack[end]))]")])
            last = pop!(tree_stack)
            push!(tree_stack, (merge(pop!(tree_stack), last)))
        end

        # Updating progress bar
        next!(prog; showvalues = [("Tree size on RAM: ", Base.format_bytes(Base.summarysize(tree_stack)))])
    end
    finish!(prog)
    # Extract the root node
    final_node = popfirst!(tree_stack)
    # JuBox.save(final_node, save_path*"leucegene_$(typeof(final_node))_samples_bt_backup.kct")
    
    # Cleanup any potential leftover leaf. Merges them linearly regardless of level.
    prog = Progress(length(tree_stack)-1, "Remaining merges to cleanup outstanding leaves"); update!(prog)
    for i in 1:length(tree_stack)-1
        pushfirst!(tree_stack, merge(popfirst!(tree_stack), popfirst!(tree_stack)))
        next!(prog)
    end
    finish!(prog)
    # JuBox.save(tree_stack[end], output)

    # Merging root node and cleanup node, which is the entirely built kct, in order of samples.
    if !isempty(tree_stack)
        final_node = merge(final_node, pop!(tree_stack))
    end

    # Saving that kct
    save(final_node, output)

end

function kct_walk(seed::String, library::String, output::String; max_recursion::Int=500, min_count::Int=2)
    kct = load(library)
    if length(kct[1][1]) > length(seed)
        println("Your seed is smaller than provided k-mers. Searching for extended seeds...")
        seeds = get_extended_seeds(seed, kct)
    else
        seeds = [seed]
    end
    fasta = ""
    for seed in seeds
        forward_root, reverse_root = build_walk_tree(seed, kct, max_recursion, min_count)
        forward_sequences = recursive_extract(forward_root.to, Vector{String}(repeat([seed], length(forward_root.to))), Dict{String, Vector{Int}}(), kct)
        reverse_sequences = recursive_extract(reverse_root.to, Vector{String}(repeat([seed], length(reverse_root.to))), Dict{String, Vector{Int}}(), kct, reverse = true)

        to_fasta(seqs, dir) = join([length(seq)>length(seed) ? ">$dir samples: $samples\n$seq\n" : "" for (seq, samples) in seqs])
        fasta *= to_fasta(forward_sequences, "forward")*to_fasta(reverse_sequences, "reverse")
    end
    open(output, "w") do file_handle
        write(file_handle, fasta)
    end
end

function build_walk_tree(seed::String, kct::Kct{N}, max_recursion::Int, min_count::Int) where {N}
    
    forward_root = Tree_Root(seed, Vector{Abstract_Node}(), [0])
    reverse_root = Tree_Root(seed, Vector{Abstract_Node}(), [0])
    forward_leaves = Vector{Abstract_Node}([forward_root])
    reverse_leaves = Vector{Abstract_Node}([reverse_root])
    max_forward_leaves = 0
    max_reverse_leaves = 0
    prog = Progress(max_recursion)
    
    for i in 1:max_recursion
        new_forward_leaves = Vector{Abstract_Node}()
        new_reverse_leaves = Vector{Abstract_Node}()
        for leaf in forward_leaves
            seq = get_seq(leaf, kct)
            alternates = [seq[2:end]*n for n in nucleotides]
            for alternate in alternates
                kct_index = findfirst(kct, DNAKmer(alternate))
                if kct_index != 0
                    samples = Vector{Int}()
                    if leaf.samples == Vector{Int}([0])
                        for (s, c) in enumerate(kct[kct_index][2])
                            if c >= min_count
                                push!(samples, s)
                            end
                        end
                        if !isempty(samples)
                            new_node = Tree_Node(kct_index, Vector{Tree_Node}(), samples)
                            push!(leaf.to, new_node)
                            push!(new_forward_leaves, new_node)
                        end
                    else
                        for s in leaf.samples
                            if kct[kct_index][2][s] >= min_count
                                push!(samples, s)
                            end
                        end
                        if !isempty(samples)
                            new_node = Tree_Node(kct_index, Vector{Tree_Node}(), samples)
                            push!(leaf.to, new_node)
                            push!(new_forward_leaves, new_node)
                        end
                    end
                end
            end
        end
        for leaf in reverse_leaves # This is disgusting, needs a rewrite
            seq = get_seq(leaf, kct)
            alternates = [n*seq[1:end-1] for n in nucleotides]
            for alternate in alternates
                kct_index = findfirst(kct, DNAKmer(alternate))
                if kct_index != 0
                    samples = Vector{Int}()
                    if leaf.samples == Vector{Int}([0])
                        for (s, c) in enumerate(kct[kct_index][2])
                            if c >= min_count
                                push!(samples, s)
                            end
                        end
                        if !isempty(samples)
                            new_node = Tree_Node(kct_index, Vector{Tree_Node}(), samples)
                            push!(leaf.to, new_node)
                            push!(new_reverse_leaves, new_node)
                        end
                    else
                        for s in leaf.samples
                            if kct[kct_index][2][s] >= min_count
                                push!(samples, s)
                            end
                        end
                        if !isempty(samples)
                            new_node = Tree_Node(kct_index, Vector{Tree_Node}(), samples)
                            push!(leaf.to, new_node)
                            push!(new_reverse_leaves, new_node)
                        end
                    end
                end
            end
        end
        if isempty(new_forward_leaves) && isempty(new_reverse_leaves) break end
        forward_leaves = new_forward_leaves
        reverse_leaves = new_reverse_leaves
        max_forward_leaves = maximum((length(forward_leaves), max_forward_leaves))
        max_reverse_leaves = maximum((length(reverse_leaves), max_reverse_leaves))
        next!(prog, showvalues = [("Current forward leaves", length(forward_leaves)),
                                  ("Current reverse leaves", length(reverse_leaves)),
                                  ("Max forward leaves", max_forward_leaves),
                                  ("Max reverse leaves", max_reverse_leaves),
                                #   ("Forward tree memory", Base.format_bytes(Base.summarysize(forward_root))),
                                #   ("Reverse tree memory", Base.format_bytes(Base.summarysize(reverse_root)))
                                  ])
    end
    finish!(prog, showvalues = [("Max forward leaves", max_forward_leaves),
                                ("Max reverse leaves", max_reverse_leaves)])
    return forward_root, reverse_root
end

function recursive_extract(current_leaves::Vector{Abstract_Node}, sequences::Vector{String}, finished::Dict{String, Vector{Int}}, kct::Kct{N}; reverse::Bool=false) where {N}
    next_leaves = Vector{Abstract_Node}()
    next_sequences = Vector{String}()
    # println(all([isempty(leaf.to) for leaf in current_leaves]))
    
    for (i, (leaf, sequence)) in enumerate(zip(current_leaves, sequences))
        # println(leaf)
        if !reverse
            to_add = get_seq(leaf, kct)[end]
            sequences[i] = sequence*to_add
        else
            to_add = get_seq(leaf, kct)[1]
            sequences[i] = to_add*sequence
        end
        if isempty(leaf.to)
            if !haskey(finished, sequences[i])
                finished[sequences[i]] = Vector{Int}()
            end
            push!(finished[sequences[i]], leaf.samples...)
        end
        # println(leaf.to)
        for next in leaf.to
            push!(next_leaves, next)
            push!(next_sequences, sequences[i])
        end
    end
    if all([isempty(leaf.to) for leaf in current_leaves]) return finished end
    return recursive_extract(next_leaves, next_sequences, finished, kct, reverse=reverse)
end

function get_extended_seeds(seed::String, kct::Kct{N}) where {N}
    prog = Progress(length(kct.table), desc="Searching for extended seeds...")
    seeds = Vector{String}()
    for k in kct.table
        if occursin(seed, String(k.seq))
            push!(seeds, k.seq)
        end
        next!(prog, showvalues=[("Extended seeds found", length(seeds))])
    end
    return seeds
end

# function prepare_leaves!(forward_root::Tree_Root, reverse_root::Tree_Root, padding::Int)
#     forward_leaves = Vector{Abstract_Node}([forward_root])
#     reverse_leaves = Vector{Abstract_Node}([reverse_root])
#     for i in 1:padding
#         new_forward_leaves =  Vector{Abstract_Node}()
#         new_reverse_leaves = Vector{Abstract_Node}()
#         for (f_leaf, r_leaf) in zip(forward_leaves, reverse_leaves)
#             for n in nucleotides
#                 new_f_leaf = Tree_Root(f_leaf.seed*n, Vector{Abstract_Node}(), 0)
#                 push!(f_leaf.to, new_f_leaf)
#                 push!(new_forward_leaves, new_f_leaf)
#                 new_r_leaf = Tree_Root(n*r_leaf.seed, Vector{Abstract_Node}(), 0)
#                 push!(r_leaf.to, new_r_leaf)
#                 push!(new_reverse_leaves, new_r_leaf)
#             end
#         end
#         forward_leaves = new_forward_leaves
#         reverse_leaves = new_reverse_leaves
#     end
#     return forward_leaves, reverse_leaves
# end

function main() 
parsed_args = parse_commands()
    if parsed_args["%COMMAND%"] == "walk"
        sub = parsed_args["walk"]
        kct_walk(sub["seed"], sub["library"], sub["output"], max_recursion = sub["max_recursion"], min_count = sub["min_count"])
    elseif parsed_args["%COMMAND%"] == "prepare"
        sub = parsed_args["prepare"]
        jf_files = readdir(sub["data"]; join = true)
        build_kct(jf_files, sub["output"], sub["jellyfish"], k=sub["k"])
    end
end

main()