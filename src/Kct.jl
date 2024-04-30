using Kmers
using JSON
using ElasticArrays
using ProgressMeter
include("parallel_sort.jl")

export Kct, Kct, load, save, merge

struct KRecord
    seq::DNA27mer
    pool_id::UInt32
end

Base.isless(a::KRecord, b::KRecord) = a.seq < b.seq
Base.isless(a::DNA27mer, b::KRecord) = a < b.seq
Base.isless(a::KRecord, b::DNA27mer) = a.seq < b
Base.:(==)(a::KRecord, b::KRecord) = a.seq == b.seq

const mask_big = 0x00000001 << 27
isbig(x::UInt32) = x & mask_big != 0
makebig(x::Integer) = UInt32(x) | mask_big
idbig(x::UInt32) = x & (~mask_big)

struct Kct{N}
    table::Vector{KRecord}
    idx::Vector{UnitRange{Int64}}
    small::ElasticMatrix{UInt8}
    big::ElasticMatrix{UInt32}
end

Kct(N, t=0, s=0, b=0) = Kct{N}(
    Vector{KRecord}(undef, t),
    fill(0:-1, 4^11),
    Matrix{UInt8}(undef, N, s),
    Matrix{UInt32}(undef, N, b)
)

function Base.push!(kct::Kct{N}, count::Matrix{UInt32}) where {N}
    if max(count...) ≤ typemax(UInt8)
        append!(kct.small, count)
        (_, id) = UInt32.(size(kct.small))
    else
        append!(kct.big, count)
        (_, id) = UInt32.(size(kct.big))
        id |= mask_big
    end
    return id
end

function Base.push!(kct::Kct{N}, count::Matrix{UInt8}) where {N}
    append!(kct.small, count)
    (_, id) = UInt32.(size(kct.small))
    return id
end

function Base.push!(kct::Kct{N}, seq::DNA27mer, count::Matrix{UInt32}) where {N}
    id = push!(kct, count)
    push!(kct.table, KRecord(seq, id))
    return id
end

function Base.getindex(kct::Kct{N}, i::Integer) where {N}
    rec = kct.table[i]
    if mask_big & rec.pool_id != 0
        tmp_t = convert(Vector{UInt32}, kct.big[:, rec.pool_id & (~mask_big)])
    else
        tmp_t = convert(Vector{UInt32}, kct.small[:, rec.pool_id])
    end
    return rec.seq, tmp_t
end


function Base.getindex(kct::Kct{N}, idx::AbstractArray{M}) where {N, M<:Integer}
    return [kct[i] for i in idx]
end

#= 
Deletes an element in the table, but does not check if the corresponding count vector
is no longer being pointed at by any element in the table. This can lead to garbage 
data left in memory until the whole Kct is freed.
Do not rely on this to make a Kct smaller. This is just a convenience thing
for when we need to work with indexes on a Kct.
=#
function Base.deleteat!(kct::Kct{N}, idx::Int64) where {N}
    deleteat!(kct.table, idx)
end

function Base.deleteat!(kct::Kct{N}, idx::AbstractArray{M}) where {N, M<:Integer}
    for i in idx
        deleteat!(kct, i)
    end
end

#= 
deleteat scales very powerly with big kcts, since it requires all elements beyond the
deleted one to be moved in memory.
quick_delete instead moves the last element to the position of the elmenent to delete
and then cuts of the tail of the table, meaning no re-assignation nescessary, which
means this is in constant time, regardless of kct size (which is pretty nice).
However, this does NOT preserve order in the table, which means this should
never be used before the kct is finished merging with others, or if you
want to be able to search k-mers using dichotomic search.
=#
function quick_delete!(kct::Kct{N}, idx::Int64) where {N}
    kct.table[idx] = kct.table[end]
    deleteat!(kct.table, length(kct.table))
end

function quick_delete!(kct::Kct{N}, idx::AbstractArray{M}) where {N, M<:Integer}
    for i in idx
        quick_delete!(kct, i)
    end
end

function Base.iterate(kct::Kct, state=1)
    state > length(kct) && return nothing
    return (kct[state], state+1)
end

Base.length(kct::Kct) = length(kct.table)

# Base.getindex(kct::Kct{N,T}, v::Vector{U}) where {N,T,U<:Integer} = [kct[i] for i in v]

# function Base.getindex(kct::Kct{N, T}, key::DNA27mer{K})::Tuple{DNA27mer{K}, NTuple{N, UInt32}} where {N, T}
#     idx_key = (key.data[1] >> 40) + 1
#     r = kct.idx[idx_key]
#     t = @view kct.table[r]
#     i = searchsortedfirst(t, key)
#     if i > length(t) || t[i].seq != key
#         return (DNA27mer{K}("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), NTuple{N, UInt32}(zeros(UInt32, N)))
#     end
#     first(kct.pool[t[i].pool_id]) == typemax(T) && return kct.next[key]
#     return (t[i].seq, convert(NTuple{N, UInt32}, kct.pool[t[i].pool_id]))
# end

# Base.getindex(kct::Kct, key::String) = kct[DNA27mer{K}(key)]

function Base.show(io::IO, t::Tuple{DNA27mer, NTuple{N, UInt32}}) where {N}
    print(io, t[1], ": ", join(string.(t[2]), ", "))

end

function Base.sort!(kct::Kct)
    psort!(kct.table)
end

function computeIndex!(kct::Kct)
    start = 1
    last_key = 0x0000000000000000
    for i=1:length(kct.table)
        key = kct.table[i].seq.data[1] >> 40
        if key > last_key
            kct.idx[last_key+1] = start:i
            start = i
            last_key = key
        end
    end
end

function _dedup(kct, nb)
    o = psortperm(nb)
    last = 0
    id = 0
    buf = Matrix{UInt32}(undef, 1, 1)
    for i in o
        if nb[i] != last
            buf[1, 1] = nb[i]
            id = push!(kct, buf)
        end
        kct.table[i] = KRecord(kct.table[i].seq, id)
        last = nb[i]
    end

    return kct
end

"""
    Kct(fn::String)

Parse a Jellyfish K-mer count table into a sorted array of DNA27mer. No
verification is made on the header, use with caution. 
"""
function Kct(fn::String; k::Int=27)
    f = open(fn, "r")
    offset = parse(Int, readuntil(f, "{"))
    json_start = position(f)
    seek(f, position(f) - 1) # Need to grab the {
    header = JSON.parse(f)
    println(header)
    seek(f, offset + json_start - 1) 

    nb = Vector{UInt32}()
    kct = Kct(1)

    while(!eof(f))
        tmp_seq = DNA27mer((read(f, UInt64),))
        tmp_nb = read(f, UInt32)
        # println(tmp_seq, " : ", tmp_nb)
        # return
        push!(kct.table, KRecord(tmp_seq, UInt32(0)))
        push!(nb, tmp_nb)
    end

    _dedup(kct, nb)
    
    sort!(kct)

    # computeIndex!(kct)

    return kct
end

# uses the dump instead of the binary. tedious and inelegant, but a good backup_kct
# in case reading the binary fails
# In our install jf executable is at /soft/bioinfo/linux_RH7/python-2.7.6/bin/jellyfish
function backup_kct(fn::String, jf::String; k::Int64=27)
    io = IOBuffer()
    cmd = pipeline(`$jf dump $fn`; stdout=io, stderr=devnull)
    run(cmd)
    data = IOBuffer(String(take!(io)))  # OUCH!!!
    nb = Vector{UInt32}()
    kct = Kct(1)
    prog = ProgressUnknown()
    while(!eof(data))
        tmp_nb = parse(UInt32, readline(data)[2:end])
        tmp_seq = DNA27mer(readline(data))
        push!(kct.table, KRecord(tmp_seq, UInt32(0)))
        push!(nb, tmp_nb)
        next!(prog)
    end

    _dedup(kct, nb)
    
    sort!(kct)

    # computeIndex!(kct)

    return kct
end

function Base.write(s::IO, kct::Kct{N}) where {N}
    write(s, N)
    write(s, length(kct.table))
    write(s, size(kct.small)[2])
    write(s, size(kct.big)[2])
    write(s, kct.table)
    write(s, kct.small)
    write(s, kct.big)
end

function Base.read(s::IO, ::Type{Kct})
    N = read(s, Int)
    l = read(s, Int)
    ls = read(s, Int)
    lb = read(s, Int)
    kct = Kct(N, l, ls, lb)
    read!(s, kct.table)
    read!(s, kct.small)
    read!(s, kct.big)
    computeIndex!(kct)
    return kct
end

function Base.findfirst(kct::Kct{N}, key::DNA27mer) where {N}
    idx_key = (key.data[1] >> 40) + 1
    r = kct.idx[idx_key]
    t = @view kct.table[r]
    i = searchsortedfirst(t, key)
    if i > length(t) || t[i].seq != key
        return 0
    end
    return r[1] + i-1
end

"""
    save(kct::Kct, fn::String)

Saves a K-mer count table into a binary format file.
"""
function save(kct::Kct{N}, fn::String) where {N}
    f = open(fn, "w")
    write(f, kct)
    close(f)
end

"""
    load(fn::String)::Kct

Loads a K-mer count table from a binary format file.
"""
function load(fn::String)::Kct
    f = open(fn, "r")
    return read(f, Kct)
end

# function _getpool(kct::Kct, pool_id::UInt32)::Union{ElasticVector{UInt8}, ElasticVector{UInt32}}
#     isbig(pool_id) && return kct.big[:,idbig(pool_id)]
#     return kct.small[:,pool_id]
# end

function _mergededup(kct::Kct, a::Kct{M}, b::Kct{N}, p) where {M,N}
    o = psortperm(p)
    last::Tuple{UInt32, UInt32} = (0, 0)
    id::UInt32 = 0
    for i in o
        if p[i] != last
            ia, ib = p[i]
            if isbig(ia) || isbig(ib)
                bbuf = zeros(UInt32, M+N, 1)
                if isbig(ia)
                    bbuf[1:M,1] = a.big[:,idbig(ia)]
                elseif ia != 0
                    bbuf[1:M,1] = a.small[:,ia]
                end
                if isbig(ib)
                    bbuf[(M+1):(M+N),1] = b.big[:,idbig(ib)]
                elseif ib != 0
                    bbuf[(M+1):(M+N),1] = b.small[:,ib]
                end
                id = push!(kct, bbuf)
            else
                buf = zeros(UInt8, M+N, 1)
                if ia != 0
                    buf[1:M,1] = a.small[:,ia]
                end
                if ib != 0
                    buf[(M+1):(M+N),1] = b.small[:,ib]
                end
                id = push!(kct, buf)
            end
        end
        kct.table[i] = KRecord(kct.table[i].seq, id)
        last = p[i]
    end
end

function Base.merge(a::Kct{M}, b::Kct{N}) where {M, N}
    i = 1
    j = 1

    kct = Kct(M+N)
    p = Vector{Tuple{UInt32, UInt32}}()

    while i ≤ length(a) || j ≤ length(b)
        if j > length(b) || (i ≤ length(a) && a.table[i].seq < b.table[j].seq)
            push!(p, (a.table[i].pool_id, 0))
            push!(kct.table, KRecord(a.table[i].seq, 0))
            i += 1
        elseif i > length(a) || (j ≤ length(b) && b.table[j].seq < a.table[i].seq)
            push!(p, (0, b.table[j].pool_id))
            push!(kct.table, KRecord(b.table[j].seq, 0))
            j += 1
        else
            push!(p, (a.table[i].pool_id, b.table[j].pool_id))
            push!(kct.table, KRecord(a.table[i].seq, 0))
            i += 1
            j += 1
        end
    end

    _mergededup(kct, a, b, p)
    # return kct, p
    return kct
end