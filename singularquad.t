local io = terralib.includec("stdio.h")
--local math = require("mathfuns")

local size_t = uint64

local math = {}
math.abs = macro(function(x) return `terralib.select(x>0, x, -x) end)
math.max = macro(function(x, y) return `terralib.select(x>y, x, y) end)
math.min = macro(function(x, y) return `terralib.select(x<y, x, y) end)

local Interval = terralib.memoize(function(T)

    --interval definition in center of gravity coordinates
    --chosen because it does not depend on orientation
    local struct interval{
        center : T
        reach : T
    }

    interval.staticmethods = {}

    interval.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or interval.staticmethods[methodname]
    end

    interval.staticmethods.new = macro(function(a, b)
        return `interval{0.5*(a+b), 0.5*math.abs(b-a)}
    end)

    --testing for equality of intervals
    interval.metamethods.__eq = terra(self : interval, other : interval) : bool
        return (self.center==other.center and self.reach==other.reach)
    end

    interval.metamethods.__ne = terra(self : interval, other : interval) : bool
        return not (self==other)
    end

    interval.metamethods.__apply = terra(self : interval, x : T) : T
        return x * self.reach + (1.0 - x) * self.center
    end

    terra interval:vol()
        return 2*self.reach
    end

    interval.staticmethods.intersect = terra(self : interval, other : interval)
        if self == other then
            return self
        else
            var distance = math.abs(self.center - other.center)
            if self.reach + other.reach >= distance then --there is an intersection
                var a = math.max(self.center-self.reach, other.center-other.reach)
                var b = math.min(self.center+self.reach, other.center+other.reach)
                return interval{0.5*(a+b), 0.5*(b-a)}
            end
        end
        return interval{0,-1} --defined as the empty interval
    end

    terra interval:isempty()
        return self.reach == -1
    end

    interval.metamethods.__and = terra(self : interval, other : interval)
        return interval.intersect(self, other)
    end

    return interval
end)

local Hypercube = terralib.memoize(function(T, N)

    local interval = Interval(T)

    --T - floating point data-type
    --N - dimension of embedding
    local struct hypercube{
        I : interval[N]
        valid : bool
    }

    terra hypercube:isempty()
        return not self.valid
    end

    hypercube.staticmethods = {}

    hypercube.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or hypercube.staticmethods[methodname]
    end

    hypercube.staticmethods.new = macro(function(...)
        local args = terralib.newlist{...}
        return `hypercube{array([args]), true}
    end)

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
        return hypercube{J, true}
    end

    hypercube.metamethods.__mul = terra(a : hypercube, b : hypercube)
        var c = hypercube.intersect(a, b)
        if c:isempty() then
            return c
        else
            for i = 0, N do
                if a:issingulardir(i) and not b:issingulardir(i) then
                    c.I[i] = b.I[i]
                elseif b:issingulardir(i) and not a:issingulardir(i) then
                    c.I[i] = a.I[i]
                elseif a:issingulardir(i) and b:issingulardir(i) then
                    c.I[i] = a.I[i]
                else
                    -- not a valid product - return with empty hypercube
                    c.valid = false
                    return c
                end
            end
            return c
        end
    end

    hypercube.metamethods.__div = terra(a : hypercube, b : hypercube)
        var c = hypercube.intersect(a, b)
        if c:isempty() then
            return c
        else
            for i = 0, N do
                if a:issingulardir(i) and not b:issingulardir(i) then
                    c.I[i] = b.I[i]
                elseif b:issingulardir(i) and not a:issingulardir(i) then
                    c.I[i] = a.I[i]
                elseif a:issingulardir(i) and b:issingulardir(i) then
                    c.I[i] = a.I[i]
                else
                    -- not a valid product - return with empty hypercube
                    c.valid = false
                    return c
                end
            end
            return c
        end
    end

    return hypercube
end)



local interval = Interval(double)
local hypercube = Hypercube(double, 2)

terra main()

    var I = interval.new(0, 1)
    var J = interval.new(2, 3)
    var Z = interval.intersect(I, J)

    --io.printf("interval(%0.2f, %0.2f)\n", I.center-I.reach, I.center+I.reach)
    --io.printf("interval(%0.2f, %0.2f)\n", J.center-J.reach, J.center+J.reach)
    --io.printf("interval(%0.2f, %0.2f)\n", Z.center-Z.reach, Z.center+Z.reach)

    var cube1 = hypercube.new(interval.new(0, 2), interval.new(0, 2))
    var cube2 = hypercube.new(interval.new(3, 4), interval.new(0, 2))
    var cube = hypercube.intersect(cube1, cube2)
    --var x = cube.I
    if cube:isempty() then
        io.printf("empty cube\n")
    end

    var a = hypercube.new(interval.new(0, 2), interval.new(0, 0))
    var b = hypercube.new(interval.new(0, 0), interval.new(0, 2))
    var c = a * b
    if c:isempty() then
        io.printf("empty cube\n")
    end
    for k = 0, 2 do
        var I = c.I[k]
        io.printf("interval(%0.2f, %0.2f)\n", I.center-I.reach, I.center+I.reach)
    end
    

end
main()