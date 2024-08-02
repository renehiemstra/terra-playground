local io = terralib.includec("stdio.h")
local fun = require("fun")

local tmath = require("mathfuns")
local svector = require("svector")
local err = require("assert")

local size_t = uint64

--return a terra tuple type of length N: {T, T, ..., T}
local ntuple = function(T, N)
    local t = terralib.newlist()
    for i = 1, N do
        t:insert(T)
    end
    return tuple(unpack(t))
end

local Point = svector.StaticVector

local Interval = terralib.memoize(function(T)

    --interval definition in center of gravity coordinates
    --chosen because it does not depend on orientation
    local struct interval{
        center : T  --exact center of the interval
        reach : T   --left and right reach
        origin : T  --an arbitrary origin
    }

    interval.staticmethods = {}

    interval.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or interval.staticmethods[methodname]
    end

    interval.staticmethods.new = macro(function(a, b)
        return `interval{0.5*(a+b), 0.5*tmath.abs(b-a), a}
    end)

    terra interval:setorigin(o : T)
        self.origin = o
    end

    terra interval:getorigin()
        return self.origin
    end

    terra interval:isempty()
        return self.reach == -1
    end

    interval.metamethods.__eq = terra(self : interval, other : interval) : bool
        return (self.center==other.center and self.reach==other.reach)
    end

    interval.metamethods.__ne = terra(self : interval, other : interval) : bool
        return not (self==other)
    end

    interval.metamethods.__apply = terra(self : interval, x : T) : T
        var A = 2*self.reach * terralib.select(self.origin <= self.center, 1., -1.)
        return A * x + self.origin
    end

    terra interval:barycentriccoord(y : T)
        var A = 2*self.reach * terralib.select(self.origin <= self.center, 1., -1.)
        return (y - self.origin) / A
    end

    terra interval:vol()
        return 2*self.reach
    end

    interval.staticmethods.intersect = terra(self : interval, other : interval)
        if self == other then
            return self
        else
            var distance = tmath.abs(self.center - other.center)
            if self.reach + other.reach >= distance then --there is an intersection
                var a = tmath.max(self.center-self.reach, other.center-other.reach)
                var b = tmath.min(self.center+self.reach, other.center+other.reach)
                return interval.new(a, b)
            end
        end
        return interval{0,-1} --defined as the empty interval
    end

    return interval
end)

local Hypercube = terralib.memoize(function(T, N, D)

    local point = Point(T,N)
    local interval = Interval(T)

    --T - floating point data-type
    --N - dimension of embedding
    local struct hypercube{
        I : interval[N]
        valid : bool
        perm : uint8[N]
    }

    terra hypercube:setperm()
        var alpha, beta = 0, N-1
        --get permutation
        var I : uint8[N]
        for i = 0, N do
            if self.I[i]:vol()==0 then
                I[i] = beta
                beta = beta - 1
            else
                I[i] = alpha
                alpha = alpha + 1
            end
        end
        --get inverse permutation
        for i = 0, N do
            self.perm[I[i]] = i
        end
    end

    terra hypercube:isempty()
        return not self.valid
    end

    hypercube.staticmethods = {}

    hypercube.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or hypercube.staticmethods[methodname]
    end

    local ctor = terra(I : interval[N])
        for i = 0, N do
            err.assert(not I[i]:isempty(), "Not a valid interval.")
        end
        var cube = hypercube{I, true}
        cube:setperm()
        return cube
    end

    hypercube.staticmethods.new = macro(terralib.memoize(function(...)
        local args = terralib.newlist{...}
        --if argument is an 'interval[N]' then return ctor
        if #args==1 then
            local v = args[1]
            local t = v.tree.type
            if t:isarray() then
                assert(t.N == N and t.type==interval)
                return `ctor([v])
            end
        end
        --else try to create 'interval[N]' and return 'ctor'
        assert(#args == N)
        return quote
            var I = arrayof(interval, [args])
        in
            ctor(I)
        end
    end))

    --check if direction 'i' is a singular direction
    terra hypercube:setorigin(o : point)
        for i = 0, N do
            self.I[i]:setorigin(o(i))
        end
    end

    terra hypercube:getorigin()
        var o : point
        for i = 0, N do
            o(i) = self.I[i]:getorigin()
        end
        return o
    end

    terra hypercube:getcenter()
        var c : point
        for i = 0, N do
            c(i) = self.I[i].center
        end
        return c
    end

    --check if direction 'i' is a singular direction
    terra hypercube:issingulardir(i : int)
        return self.I[i]:vol()==0
    end

    --dimension is computed at runtime as the dimension of the embedding
    --minus the singular dimensions
    terra hypercube:dim()
        var d = N
        for i = 0, N do
            if self:issingulardir(i) then
                d = d - 1
            end
        end
        return d
    end

    --compute the length, area, volume, etc. 'vol' of a point
    --is defined as '1'
    terra hypercube:vol()
        var v = 1.0
        for i = 0, N do
            if not self:issingulardir(i) then
                v = v * self.I[i]:vol()
            end
        end
        return v
    end

    terra hypercube:barycentriccoord(p : point)
        var x : ntuple(T, N)
        var ptr = [&T](&x)
        var I : interval
        var k : uint8
        var D = self:dim()
        for i = 0, D do
            k = self.perm[i]
            I = self.I[k]
            ptr[i] = I:barycentriccoord(p(k))
        end
        for i = D, N do
            ptr[i] = 0
        end
        return x
    end

    local terra eval(self : &hypercube, x : &T)
        var d = self:dim()
        var p : point
        for i = 0, d do
            var k = self.perm[i]
            p(k) = self.I[k](x[i])
        end
        for i = d, N do
            var k = self.perm[i]
            p(k) = self.I[k].center
        end
        return p
    end

    hypercube.metamethods.__apply = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        local D = #args
        --try a 'point' or a 'tuple'
        if D==1 then
            local p = args[1]
            local typ = p.tree.type
            --try 'point'
            if typ==point then
                return `eval(&self, [&T](&[ p ]))
            else 
                --try n-tuple
                if typ.convertible == "tuple" then
                    local n = #typ.entries
                    local x = symbol(T[n])
                    local fill_array = terralib.newlist{}
                    for k = 0, n-1 do
                        fill_array:insert(quote [x][k] = p.["_"..tostring(k)] end)
                    end
                    return quote
                        err.assert(self:dim() == [n])
                        var [ x ]
                        [fill_array]
                    in
                        eval(&self, [&T](&x))
                    end
                end
            end
        end
        --try list of arguments
        local x = symbol(T[D])
        local fill_array = terralib.newlist{}
        for k,v in ipairs(args) do
            assert(v.tree.type:isprimitive() or "Not a primitive type")
            fill_array:insert(quote [x][k-1] = [ v ] end)
        end
        return quote
            var [ x ]
            err.assert(self:dim()==D)
            [ fill_array ]
        in
            eval(&self, [&T](&x))
        end
    end))

    --testing for equality of intervals
    hypercube.metamethods.__eq = terra(self : hypercube, other : hypercube) : bool
        for k = 0, N do
            if self.I[k] ~= other.I[k] then
                return false
            end
        end
        return true
    end

    hypercube.metamethods.__ne = terra(self : hypercube, other : hypercube) : bool
        return not (self==other)
    end

    hypercube.staticmethods.intersect = terra(self : hypercube, other : hypercube)
        var J : interval[N]
        for k = 0, N do
            var Z = interval.intersect(self.I[k], other.I[k])
            if Z:isempty() then
                return hypercube{J, false} --return empty (invalid) hypercube
            end
            J[k] = Z
        end
        return hypercube.new(J)
    end

    hypercube.metamethods.__mul = terra(a : hypercube, b : hypercube)
        var cube = hypercube.intersect(a, b)
        if not cube:isempty() then
            for i = 0, N do
                if a:issingulardir(i) and not b:issingulardir(i) then
                    cube.I[i] = b.I[i]
                elseif b:issingulardir(i) and not a:issingulardir(i) then
                    cube.I[i] = a.I[i]
                elseif a:issingulardir(i) and b:issingulardir(i) then
                    cube.I[i] = a.I[i]
                else
                    -- not a valid product - return with empty hypercube
                    cube.valid = false
                    return cube
                end
            end
        end
        cube:setperm()
        return cube
    end

    hypercube.metamethods.__div = terra(a : hypercube, b : hypercube)
        var cube = hypercube.intersect(a, b)
        var p = cube:getorigin()
        if cube == b then
            for i = 0, N do
                if b:issingulardir(i) then
                    cube.I[i] = a.I[i]
                else
                    cube.I[i].reach = 0
                    cube.I[i].center = cube.I[i].origin
                end
            end
        else
            cube.valid = false
        end
        cube:setperm()
        return cube
    end

    return hypercube
end)


return {
    Point = Point,
    Interval = Interval,
    Hypercube = Hypercube
}