local io = terralib.includec("stdio.h")
local err = require("assert")
local tup = require("tuple")
local tmath = require("mathfuns")
local range = require("range")
local luafun = require("fun")

local size_t = uint64

local linrange = {}
linrange.__index = linrange
linrange.__metatable = linrange

function linrange.new(n)
    local t = {}
    t.n = n
    return setmetatable(t, linrange)
end

function linrange:next(i)
    local i = i+1
    if i < self.n then
        return i
    end
end

function linrange:iterate()
    return linrange.next, self, -1
end

local function next(sizes, iter, i, k)
    i[k] = iter[k]:next(i[k])
    if i[k] then
        return true
    else
        if k>1 then
            iter[k] = linrange.new(sizes[k])
            i[k] = iter[k]:next(-1)
            return next(sizes, iter, i, k-1)
        end
    end
    return false
end

local productiter = function(...)
    local sizes = terralib.newlist{...}
    local m = #sizes
    local iter = sizes:map(function(s) return linrange.new(s) end)
    local i = sizes:map(function(s) return 0 end)
    i[m] = -1
    return function()
        if next(sizes, iter, i, m) then
            return i
        end
    end
end

local getarrayentry = function(args, multiindex)
    local expr = args
    for i, index in ipairs(multiindex) do
        local s = "_" .. index
        expr = `(expr).[s]
    end
    return expr
end

local function defaultperm(dimension)
    return terralib.newlist(
        luafun.totable(
            luafun.range(dimension, 1, -1)
        )
    )
end

local function checkperm(perm)
    assert(terralib.israwlist(perm), "ArgumentError: input should be a raw list.")
    local linrange = perm:mapi(function(i,v) return i end)
    for i,v in ipairs(perm) do
        linrange[v] = nil
    end
    assert(#linrange == 0, "ArgumentError: input list is not a valid permutation.")
end


--ArrayBase collects some common utilities used to built 
--static and dynamic Array types. The following traits
--and methods are assumed
--  Array.eltype
--  Array.ndims
--  Array.perm
--  Array.methods.length
--  Array.methods.size
local ArrayBase = function(Array)

    local T = Array.eltype
    local N = Array.ndims

    local Unitrange = range.Unitrange(size_t)

    --local terra array holding the permutation array
    local __perm = terralib.constant(terralib.new(size_t[N], Array.perm))
    local rowmajorperm = defaultperm(N) --the standard permutation
    
    --local cartesian indices
    local Indices = rowmajorperm:map(function(v) return symbol(size_t) end)

    local terra boundscheck_linear(self : &Array, index : size_t)
        err.assert(index < self:length(), "BoundsError: out of bounds.")
    end

    local terra boundscheck_cartesian(self : &Array, [Indices])
        escape
            for d = 1, N do
                local index = Indices[d]
                local message = "BoundsError: array dimension " .. tostring(d) .. " out of bounds."
                emit quote
                    err.assert([index] < self:size([d-1]), message)
                end
            end
        end
    end

    local boundscheck = macro(function(self, ...)
        if __boundscheck__ then
            local indices = terralib.newlist{...}
            if #indices == 1 then
                return `boundscheck_linear(self, [indices])
            else
                return `boundscheck_cartesian(self, [indices])
            end
        end
    end)

    local function processarraydim(array, size, k)
        if k < N then
            size[k+1] = #array[1]
            for i,v in ipairs(array) do
                assert(type(v) == "table", "ArgumentError: expected array input of dimension " .. tostring(N) ..".")
                array[i] = terralib.newlist(v)
                assert(#array[i] == size[k+1], "ArgumentError: array input size inconsistent with array dimensions.")
                processarraydim(array[i], size, k+1)
            end
        end
    end

    local function processarrayinput(array)
        array = terralib.newlist(array)
        local size = terralib.newlist()
        size[1] = #array
        processarraydim(array, size, 1)
        return size
    end


    --return linear indices product range
    local function getunitranges(K)
        return terra(self : &Array)
            var uranges : tup.ntuple(Unitrange, K)
            escape
                for k = 1, K do
                    local s = "_" .. tostring(k-1)
                    emit quote uranges.[s] = Unitrange{0, self:size([k-1])} end
                end
            end
            return uranges
        end
    end

    --return linear indices product range
    local terra getunitrange(self : &Array, i : size_t)
        return Unitrange{0, self:size(i)}
    end

    local printarray
    if N == 1 then
        printarray = function(self, name)
            return quote
                io.printf("%s = \n", [ name ])
                for i in getunitrange(&self, 0) do
                    var value = self:get(i)
                    io.printf("[%s]\n", tmath.numtostr(value))
                end
                io.printf("\n")
            end
        end
    elseif N == 2 then
        printarray = function(self, name)
            return quote 
                io.printf("%s = \n", [ name ])
                for i in getunitrange(&self, 0) do
                    io.printf("\t[")
                    for j in getunitrange(&self, 1) do
                        var value = self:get(i, j)
                        io.printf("%s\t", tmath.numtostr(value))
                    end
                    io.printf("]\n")
                end
                io.printf("\n")
            end
        end
    else
        printarray = function(self, name)
            local unitranges = getunitranges(N-2) --terra function that returns the first N-2 unitranges
            local p = rowmajorperm:filteri(function(i,v) return i <= N - 2 end)
            local ntimes = p:mapi(function(i,v) return "%d" end)
            local slice = name .."[" .. table.concat(ntimes,",") .. ", :, :] = \n"
            return quote
                var K = unitranges(&self)
                for k in range.product(unpacktuple(K)) do
                    io.printf([ slice ], unpacktuple(k))
                    for i in getunitrange(&self, N-2) do
                        io.printf("\t[")
                        for j in getunitrange(&self, N-1) do
                            var value = self:get(unpacktuple(k), i, j)
                            io.printf("%s\t", tmath.numtostr(value))
                        end
                        io.printf("]\n")
                    end
                    io.printf("\n")
                end
            end
        end
    end

    local print = macro(function(self)
        local name = self.tree.name
        return printarray(self, name)
    end)

    --element size as a terra method
    local getperm = terralib.overloadedfunction("perm", {
        terra(self : &Array, i : size_t)
            return __perm[i]
        end,
        terra(self : &Array)
            return __perm
        end
    })

    return {
        rowmajorperm = rowmajorperm,
        Indices = Indices,
        boundscheck_linear = boundscheck_linear,
        boundscheck_cartesian = boundscheck_cartesian,
        boundscheck = boundscheck,
        processarrayinput = processarrayinput,
        getperm = getperm,
        getunitranges = getunitranges,
        print = print
    }

end


return {
    getarrayentry = getarrayentry,
    productiter = productiter,
    checkperm = checkperm,
    defaultperm = defaultperm,
    ArrayBase = ArrayBase
}