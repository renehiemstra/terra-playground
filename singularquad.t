local io = terralib.includec("stdio.h")

--local math = require("mathfuns")
local err = require("assert")

local size_t = uint64

local math = {}
math.abs = macro(function(x) return `terralib.select(x>0, x, -x) end)
math.max = macro(function(x, y) return `terralib.select(x>y, x, y) end)
math.min = macro(function(x, y) return `terralib.select(x<y, x, y) end)

local Point = terralib.memoize(function(T, N)

    local struct point{
        position : T[N]
    }

    point.staticmethods = {}

    point.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or point.staticmethods[methodname]
    end

    point.staticmethods.new = macro(function(...)
        local args = terralib.newlist{...}
        assert(#args == N)
        return `point{arrayof(T,[args])}
    end)

    point.metamethods.__apply = macro(function(self, i)
        return quote
            err.assert(i >= 0 and i < N)
        in
            self.position[i]
        end
    end)

    return point
end)


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
        return x * (self.center + self.reach) + (1.-x) * (self.center - self.reach)
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
        assert(#args == N)
        return quote
            var A = hypercube{array([args]), true}
            for i = 0, N do
                if A.I[i]:isempty() then
                    A.valid = false
                    break
                end
            end
        in
            A
        end
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
        if not c:isempty() then
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
        end
        return c
    end

    hypercube.metamethods.__div = terra(a : hypercube, b : hypercube)
        var c = hypercube.intersect(a, b)
        if c == b then
            for i = 0, N do
                if b:issingulardir(i) then
                    c.I[i] = a.I[i]
                else
                    c.I[i].reach = 0
                end
            end
        else
            c.valid = false
        end
        return c
    end

    terra hypercube:print(name : rawstring)
        var I = self.I
        io.printf("%s : (%0.2f, %0.2f) - (%0.2f, %0.2f)\n", name, I[0](0), I[1](0), I[0](1), I[1](1))
    end

    return hypercube
end)

local Pyramid = terralib.memoize(function(T, N)

    local interval = Interval(T)

    --T - floating point data-type
    --N - dimension of embedding
    local struct pyramid{
        I : interval[N]
        valid : bool
    }

    pyramid.staticmethods = {}

    pyramid.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or pyramid.staticmethods[methodname]
    end

    pyramid.staticmethods.new = macro(function(...)
        local args = terralib.newlist{...}
        assert(#args == N)
    end)

end)


import "terratest/terratest"


testenv "point" do

    local point = Point(double, 3)
    
    testset "new, apply" do
        terracode
            var p = point.new(1, 3, 4)
            --p(2) = 5.0
        end
        test p(0)==1
        test p(1)==3
        test p(2)==5
    end

end

testenv "interval" do

    local interval = Interval(double)
    
    testset "new, vol, apply" do
        terracode
            var I = interval.new(1, 3)
        end
        test I.center==2
        test I.reach==1
        test I:vol()==2
        test I(0)==1
        test I(1)==3
        test I(0.5)==I.center
    end

    testset "empty intersection" do
        terracode
            var I = interval.new(0, 1)
            var J = interval.new(2, 3)
            var Z = interval.intersect(I, J)
        end
        test Z.reach == -1
        test Z:isempty()
    end

    testset "nonempty intersection - interval" do
        terracode
            var I = interval.new(0, 2)
            var J = interval.new(1, 3)
            var Z = interval.intersect(I, J)
        end
        test Z==interval.new(1, 2)
        test Z.center==1.5
        test Z.reach==0.5
    end

    testset "nonempty intersection - point" do
        terracode
            var I = interval.new(0, 1)
            var J = interval.new(1, 2)
            var Z = interval.intersect(I, J)
        end
        test Z==interval.new(1, 1)
        test Z.center==1
        test Z.reach==0
    end

    testset "comparison" do
        terracode
            var I = interval.new(0, 2)
            var J = interval.new(0, 2)
            var Z = interval.new(1, 2)
        end
        test I == J
        test I ~= Z
    end
end


testenv "hypercube" do

    local interval = Interval(double)
    local hypercube = Hypercube(double, 2)
    
    testset "new" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(0, 0))
            var C = hypercube.new(interval.new(0, 0), interval.new(0, 0))
        end
        test A:isempty()==false and B:isempty()==false and C:isempty()==false
        test A:dim()==2 and A:issingulardir(0)==false and A:issingulardir(1)==false
        test B:dim()==1 and B:issingulardir(0)==false and B:issingulardir(1)==true
        test C:dim()==0 and C:issingulardir(0)==true  and C:issingulardir(1)==true
    end

    testset "empty intersection" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(3, 4), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()
    end

    testset "nonempty intersection - surface" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(1, 3), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(1, 2), interval.new(0, 2))
    end

    testset "nonempty intersection - line" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(2, 3), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(2, 2), interval.new(0, 2))
    end

    testset "nonempty intersection - point" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(2, 3), interval.new(2, 3))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(2, 2), interval.new(2, 2))
    end

    testset "product" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 0))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 2))
            var C = A * B
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(0, 2), interval.new(0, 2))
    end

    testset "division" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 0))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 2))
            var C = A * B
        end
        test (C / A) * A == C
        test (C / B) * B == C
    end
end

local interval = Interval(double)
local hypercube = Hypercube(double, 2)

terra main()
    var A = hypercube.new(interval.new(0, 2), interval.new(0, 0))
    var B = hypercube.new(interval.new(0, 0), interval.new(0, 2))
    var C = A * B
    var D = C / A
    
    A:print("A")
    B:print("B")
    C:print("C")
    D:print("D")
end
main()