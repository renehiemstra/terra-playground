local io = terralib.includec("stdio.h")

local math = require("mathfuns")
local err = require("assert")

local size_t = uint64

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
        return `self.position[i]
    end)

    point.metamethods.__eq = terra(self : point, other : point) : bool
        for i = 0, N do
            if self.position[i] ~= other.position[i] then
                return false
            end
        end
        return true
    end

    point.metamethods.__ne = terra(self : point, other : point) : bool
        return not (self==other)
    end

    return point
end)


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
        return `interval{0.5*(a+b), 0.5*math.abs(b-a), 0.5*(a+b)}
    end)

    terra interval:setorigin(o : T)
        self.origin = o
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
        for i = 0, N do
            if self.I[i]:vol()==0 then
                self.perm[i] = beta
                beta = beta - 1
            else
                self.perm[i] = alpha
                alpha = alpha + 1
            end
        end
    end

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
            if not A:isempty() then
                A:setperm()
            end
        in
            A
        end
    end)

    --check if direction 'i' is a singular direction
    terra hypercube:setorigin(o : point)
        for i = 0, N do
            self.I[i]:setorigin(o.position[i])
        end
    end

    terra hypercube:getcenter()
        var c : point
        for i = 0, N do
            c.position[i] = self.I[i].center
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

    hypercube.metamethods.__apply = macro(function(self, ...)
        local args = terralib.newlist{...}
        local D = #args
        local eval = terralib.newlist()
        local I = symbol(interval)
        local p = symbol(point)
        local k = symbol(int)
        for i = 0, D-1 do
            local x = args[i+1]
            eval:insert(quote
                [k] = self.perm[i]
                [I] = self.I[k]
                [p].position[k] = I([x])
            end)
        end
        for i = D, N-1 do
            eval:insert(quote
                [k] = self.perm[i]
                [I] = self.I[k]
                [p].position[k] = I.center
            end)
        end
        return quote
            var [ k ] = 0
            var [ I ] = interval{}
            var [ p ] = point{}
            err.assert(D == self:dim())
            [eval]
        in
            p
        end
    end)

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

    local point = Point(T,N)
    local hypercube = Hypercube(T, N)
    --dimension 'D' of base should be larger than
    --zero and smaller than embedding dimension N
    --but D is a runtime parameter

    --T - floating point data-type
    --N - dimension of embedding
    local struct pyramid{
        base : hypercube
        apex : point
    }

    pyramid.staticmethods = {}

    pyramid.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or pyramid.staticmethods[methodname]
    end

    pyramid.staticmethods.new = terra(base : hypercube, apex : point)
        var D = base:dim()
        err.assert(D > 0 and D < N)
        return pyramid{base, apex}
    end

    terra pyramid:dim()
        return self.base.dim()+1
    end

    terra pyramid:height()
        return 1.0
    end

    terra pyramid:vol()
        return self:height() * self.base:vol() / self:dim()
    end

    pyramid.metamethods.__apply = terra(self : &pyramid)
    end

end)


import "terratest/terratest"


testenv "point" do

    local point = Point(double, 3)
    
    testset "new, apply" do
        terracode
            var p = point.new(1, 3, 4)
        end
        test p(0)==1
        test p(1)==3
        test p(2)==4
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
        test I.origin == I.center
        test I:vol()==2
    end

    testset "evaluation" do
        terracode
            var I = interval.new(1, 3)
            I.origin = 1.2
            var J = interval.new(1, 3)
            J.origin = 2.8
        end
        --test I(-0.1)==1 and I(0.9)==3
        test math.isapprox(I(-0.1), 1.0, 1e-15) 
        test math.isapprox(I(0.9), 3.0, 1e-15)
        test math.isapprox(J(-0.1), 3.0, 1e-15) 
        test math.isapprox(J(0.9), 1.0, 1e-15)
    end

    testset "barycentric coordinate - matching orientation" do
        terracode
            var I = interval.new(1, 3)
            I:setorigin(1)
        end
        test I:barycentriccoord(1)==0
        test I:barycentriccoord(2)==0.5
        test I:barycentriccoord(3)==1
    end

    testset "barycentric coordinate - reverse orientation" do
        terracode
            var I = interval.new(1, 3)
            I:setorigin(3)
        end
        test I:barycentriccoord(3)==0
        test I:barycentriccoord(2)==0.5
        test I:barycentriccoord(1)==1
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


testenv "hypercube - 2" do

    local point = Point(double, 2)
    local interval = Interval(double)
    local hypercube = Hypercube(double, 2)
    
    testset "new, dim, vol" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(0, 0))
            var C = hypercube.new(interval.new(0, 0), interval.new(0, 0))
        end
        test A:isempty()==false and B:isempty()==false and C:isempty()==false
        test A:dim()==2 and A:issingulardir(0)==false and A:issingulardir(1)==false
        test B:dim()==1 and B:issingulardir(0)==false and B:issingulardir(1)==true
        test C:dim()==0 and C:issingulardir(0)==true  and C:issingulardir(1)==true
        test A:vol()==4 and B:vol()==2 and C:vol()==1
    end

    testset "permutation" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5))
            var B = hypercube.new(interval.new(0, 2), interval.new(1, 1))
            var C = hypercube.new(interval.new(0, 0), interval.new(1, 2))
        end
        test A.perm[0] == 0 and A.perm[1] == 1
        test B.perm[0] == 0 and B.perm[1] == 1
        test C.perm[0] == 1 and C.perm[1] == 0
    end

    testset "square evaluation - origin==center" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5))
        end
        test A(0,0) == A:getcenter()
        test A(0.5,0.5) == point.new(2,5)
        test A(-0.5,-0.5) == point.new(0,1)
    end

    testset "square evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5))
            A:setorigin(point.new(2,5))
        end
        test A(0,0) == point.new(2,5)
        test A(1,1) == point.new(0,1)
    end

    testset "segment evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(1, 3), interval.new(1, 1))
            A:setorigin(point.new(3,1))
        end
        test A(0) == point.new(3,1)
        test A(1) == point.new(1,1)
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

testenv "hypercube - 3" do

    local point = Point(double, 3)
    local interval = Interval(double)
    local hypercube = Hypercube(double, 3)
    
    testset "new, dim, vol" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 0))
            var C = hypercube.new(interval.new(0, 2), interval.new(0, 0), interval.new(0, 0))
            var D = hypercube.new(interval.new(0, 0), interval.new(0, 0), interval.new(0, 0))
        end
        test A:isempty()==false and B:isempty()==false and C:isempty()==false and D:isempty()==false
        test A:dim()==3 and A:issingulardir(0)==false and A:issingulardir(1)==false and A:issingulardir(2)==false
        test B:dim()==2 and B:issingulardir(0)==false and B:issingulardir(1)==false and B:issingulardir(2)==true
        test C:dim()==1 and C:issingulardir(0)==false and C:issingulardir(1)==true  and C:issingulardir(2)==true
        test D:dim()==0 and D:issingulardir(0)==true  and D:issingulardir(1)==true  and D:issingulardir(2)==true
        test A:vol()==8 and B:vol()==4 and C:vol()==2 and D:vol()==1
    end

    testset "permutation" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 2), interval.new(0, 2))
            var C = hypercube.new(interval.new(0, 0), interval.new(0, 1), interval.new(0, 0))
            var D = hypercube.new(interval.new(0, 0), interval.new(0, 0), interval.new(0, 0))
        end
        test A.perm[0] == 0 and A.perm[1] == 1 and A.perm[2] == 2
        test B.perm[0] == 2 and B.perm[1] == 0 and B.perm[2] == 1
        test C.perm[0] == 2 and C.perm[1] == 0 and C.perm[2] == 1
        test D.perm[0] == 2 and D.perm[1] == 1 and D.perm[2] == 0
    end

    testset "cube evaluation - origin==center" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5), interval.new(2, 4))
        end
        test A(0,0,0) == A:getcenter()
        test A(0.5,0.5,0.5) == point.new(2,5,4)
        test A(-0.5,-0.5,-0.5) == point.new(0,1,2)
    end

    testset "cube evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5), interval.new(2, 4))
            A:setorigin(point.new(2,5,4))
        end
        test A(0,0,0) == point.new(2,5,4)
        test A(1,1,1) == point.new(0,1,2)
    end

    testset "surface evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 1), interval.new(2, 4))
            A:setorigin(point.new(0,1,2))
        end
        test A(0,0) == point.new(0,1,2)
        test A(1,0) == point.new(2,1,2)
    end

    testset "segment evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(2, 2), interval.new(1, 1), interval.new(1, 3))
            A:setorigin(point.new(2,1,3))
        end
        test A(0) == point.new(2,1,3)
        test A(1) == point.new(2,1,1)
    end

    testset "empty intersection" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(3, 4), interval.new(0, 2), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()
    end

    testset "nonempty intersection - cube" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(1, 3), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(0, 2), interval.new(1, 2), interval.new(0, 2))
    end

    testset "nonempty intersection - surface" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(2, 3))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(2, 2))
    end

    testset "nonempty intersection - line" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(2, 3), interval.new(0, 2), interval.new(2, 3))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(2, 2), interval.new(0, 2), interval.new(2, 2))
    end

    testset "nonempty intersection - point" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(2, 3), interval.new(2, 3), interval.new(2, 3))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(2, 2), interval.new(2, 2), interval.new(2, 2))
    end

    testset "product" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 0))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 0), interval.new(0, 1))
            var C = A * B
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 1))
    end

    testset "division" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 0))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 0), interval.new(0, 1))
            var C = A * B
        end
        test (C / A) * A == C
        test (C / B) * B == C
    end
end

local interval = Interval(double)
local hypercube = Hypercube(double, 2)
local point = Point(double, 2)

terra main()
    
    var A = hypercube.new(interval.new(0, 2), interval.new(0, 0))
    A:setorigin(point.new(2,0))
    io.printf("evaluate A(x, y) = (%0.2f, %0.2f)\n", A(0))
end
main()