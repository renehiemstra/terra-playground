local io = terralib.includec("stdio.h")
local err = require("assert")
local tmath = require("mathfuns")

local size_t = uint64


local StackBase = terralib.memoize(function(S, T, N)

    S.eltype = T

    S.methods.size = terra(self : &S) : size_t
        return N
    end

    S.methods.get = terra(self : &S, i : size_t) : T
        err.assert(i < N)
        return self.data[i]
    end

    S.methods.set = terra(self : &S, i : size_t, value : T)
        err.assert(i < N)
        self.data[i] = value
    end

    S.metamethods.__apply = macro(function(self, i)
        return quote
            err.assert(i < N)
        in
            self.data[i]
        end
    end)

    return S
end)

local VectorBase = terralib.memoize(function(V, T, N)

    --add all functionality of StackBase
    StackBase(V,T,N)

    V.staticmethods = {}

    V.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or V.staticmethods[methodname]
    end

    V.staticmethods.fill = terra(value : T)
        var v : V
        v.simd = value
        return v
    end

    V.staticmethods.zeros = terra()
        return V.fill(T(0))
    end

    V.staticmethods.ones = terra()
        return V.fill(T(1))
    end

    V.staticmethods.new = terra()
        return V.fill(T(0))
    end

    V.staticmethods.from = macro(function(...)
        local args = terralib.newlist{...}
        assert(#args == N, "Length of input list does not match static dimension")
        local eval = terralib.newlist{}
        local x = symbol(V)
        for k, v in ipairs(args) do
            eval:insert(quote [x].data[ [k-1] ] = [ v ] end)
        end
        return quote
            var [x]
            [eval]
        in
            x
        end
    end)

    V.metamethods.__eq = terra(self : V, other : V)
        var v = self.simd == other.simd
        for i = 0, N do
            if not v[i] then return false end
        end
        return true
    end

    V.metamethods.__ne = terra(self : V, other : V)
        return (self==other)==false
    end

    local function getunderlying(variable)
        local typ = variable.tree.type
        if typ:isprimitive() then
            return `variable
        elseif typ==V then
            return `variable.simd
        end
    end

    V.metamethods.__add = macro(function(self, other)
        local v1, v2 = getunderlying(self), getunderlying(other)
        return `V{[v1] + [v2]}
    end)

    V.metamethods.__sub = macro(function(self, other)
        local v1, v2 = getunderlying(self), getunderlying(other)
        return `V{v1 - v2}
    end)
    
    V.metamethods.__mul = macro(function(self, other)
        local v1, v2 = getunderlying(self), getunderlying(other)
        return `V{[v1] * [v2]}
    end)

    --doesn't work on integral types (at least on mac)
    if not T:isintegral() then
        V.metamethods.__div = macro(function(self, other)
            local v1, v2 = getunderlying(self), getunderlying(other)
            return `V{[v1] / [v2]}
        end)
    end

    local __sum = macro(function(self)
        local eval = terralib.newlist()
        local s = symbol(T)
        for i = 0, N-1 do
            eval:insert(quote [ s ] = [ s ] + self.data[i] end)
        end
        return quote
            var [s] = 0
            [eval]
        in
            [s]
        end
    end)

    terra V:sum()
        return __sum(self)
    end

    terra V:norm2()
        return __sum(@self * @self)
    end

    terra V:norm()
        return tmath.sqrt([double](self:norm2()))
    end
    
    return V
end)


local StaticVector = terralib.memoize(function(T, N)

    --requesting less than 32 bytes results in a bug on macos
    --so for now we specifically request 32 bytes or more.
    local M = N
    if (sizeof(T)*N < 32) then M = (32 / sizeof(T)) end
    local SIMD = vector(T,M)

    --for the union make sure M corresponds to sizeof(SIMD) / sizeof(T) 
    --because sizeof(SIMD) is always a multiple of 2
    M = sizeof(SIMD) / sizeof(T)

    local Base = function(V) VectorBase(V,T,N) end

    local struct svector(Base){
        union{
            simd : SIMD
            data : T[M]
        }
    }

    svector.methods.print = terra(self : &svector, name : rawstring)
        io.printf("%s = (", name)
        for i = 0, N-1 do
            io.printf("%0.2f, ", self(i))
        end
        io.printf("%0.2f)\n", self(N-1))
    end
    
    return svector
end)


return {
    StaticVector = StaticVector
}