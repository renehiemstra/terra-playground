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

    interval.metamethods.__cast = terralib.memoize(function(from, to, exp)
        if from==interval and to:isprimitive() then
            local S = to
            return quote
                var v = exp
                err.assert(v.reach==0, "Non-singular interval cannot be cast to a real number")
            in
                [S](v.center)
            end
        end
        if from:isprimitive() and to==interval then
            return quote
                var v = exp
            in
                interval.new(v, v)
            end
        end
    end)

    return interval
end)

local Cube = terralib.memoize(function(T, N, D)
    assert(D <= N and D >= 0)
    local point = Point(T,N)
    local point_d = Point(T,D)
    local interval = Interval(T)

    --T - floating point data-type
    --N - dimension of embedding
    local struct cube{
        I : interval[D]
        J : T[N-D]
        perm : uint8[N]
        valid : bool
    }

    terra cube:isempty()
        return not self.valid
    end

    cube.staticmethods = {}

    cube.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or cube.staticmethods[methodname]
    end

    cube.staticmethods.new = macro(terralib.memoize(function(...)
        local args = terralib.newlist{...}
        assert(#args == N and "Number of arguments inconsistent with dimension.")
        local intervals = terralib.newlist()
        local positions = terralib.newlist()
        return quote
            var invperm : uint8[N]
            var perm : uint8[N]
            var I : interval[D]
            var J : T[N-D]
            escape 
                --compute inverse permutation
                local alpha, beta = 0, N
                for i, v in ipairs(args) do
                    if v.tree.type==interval then
                        emit quote invperm[i-1] = [alpha] end
                        alpha = alpha + 1
                        intervals:insert(v)
                    elseif v.tree.type:isprimitive() then
                        beta = beta - 1
                        emit quote invperm[i-1] = [beta] end
                        positions:insert(v)
                    end
                end
                --check argument counts
                assert(#intervals == D and #positions == N-D and "number of intervals inconsistent with dimension")
                --compute permutation
                for i = 0, N-1 do
                    emit quote perm[ invperm[i] ] = i end
                end
                --fill intervals array
                for i, v in ipairs(intervals) do
                    emit quote I[i-1] = v end
                end
                --fill positions array
                for i, v in ipairs(positions) do
                    emit quote J[i-1] = v end
                end
            end
        in
            cube{I, J, perm, true}
        end
    end))

    terra cube:setorigin(o : point)
        var k = 0
        escape
            for i = 0, D-1 do
                emit quote 
                    k = self.perm[i]
                    self.I[i]:setorigin(o(k))
                end
            end
        end
    end

    --check if direction 'i' is a singular direction
    terra cube:issingulardir(i : int)
        return self.perm[i] >= D
    end

    terra cube:getinterval(i : int) : interval
        var k = self.perm[i]
        if k < D then
            return self.I[k]
        else
            return self.J[k-D]
        end
    end
    
    terra cube:dim() return D end

    --compute the length, area, volume, etc. 'vol' of a point
    --is defined as '1'
    terra cube:vol()
        var v = 1.0
        escape
            for i = 0, D-1 do
                emit quote v = v * self.I[i]:vol() end
            end
        end
        return v
    end

    local terra eval(self : &cube, x : &point_d)
        var y : point
        var k : uint8
        escape
            for i = 0, D-1 do
                emit quote
                    k = self.perm[i]
                    y(k) = self.I[i](x(i))
                end
            end
            for i = D, N-1 do
                emit quote
                    k = self.perm[i]
                    y(k) = self.J[i-D]
                end
            end
        end
        return y
    end

    cube.metamethods.__apply = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        --first try a point or a tuple
        if #args==1 and args[1].tree.type==point_d then
            local v = args[1]
            return `eval(&self, &[v])
        else
            return quote
                var x = point_d.from([args])
            in
                eval(&self, &x)
            end
        end
    end))

    terra cube:barycentriccoord(p : point)
        var x : point_d
        var k : uint8
        escape
            for i = 0, D-1 do
                emit quote
                    k = self.perm[i]
                    x(i) = self.I[k]:barycentriccoord(p(k))
                end
            end
        end
        return x
    end
       
    return cube
end)


local Hypercube = terralib.memoize(function(T, N)

    --get all dimenions of Cube
    local cube = terralib.newlist{}
    for D = 0, N do
        cube[D] = Cube(T,N,D)
    end

    --attach all methods to this table
    local hypercube = {}

    --add methods for overloaded function 'intersection'
    hypercube.intersect = macro(function(a, b)

    end)

end)




--return {
--    Point = Point,
--    Interval = Interval,
--    Cube = Cube
--}

import "terratest/terratest"


testenv "interval" do

    local interval = Interval(double)
    
    testset "new, vol, apply" do
        terracode
            var I = interval.new(1, 3)
        end
        test I.center==2
        test I.reach==1
        test I.origin == 1
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
        test tmath.isapprox(I(-0.1), 1.0, 1e-15) 
        test tmath.isapprox(I(0.9), 3.0, 1e-15)
        test tmath.isapprox(J(-0.1), 3.0, 1e-15) 
        test tmath.isapprox(J(0.9), 1.0, 1e-15)
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

    testset "casting" do
        terracode
            var I = interval.new(1, 1)
        end
        test [double](I)==1
        test [interval](1)==I
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

testenv "cube - N = 2" do

    local T = double
    local point_1d = Point(T, 1)
    local point_2d = Point(T, 2)
    local interval = Interval(T)
    local square = Cube(T, 2, 2)
    local line = Cube(T, 2, 1)

    terracode
        var A = square.new(interval.new(0, 2), interval.new(1, 5))
        var B = line.new(interval.new(0, 2), 1)
        var C = line.new(0, interval.new(1, 2))
    end
    
    testset "new, dim, vol" do
        test A:isempty()==false and B:isempty()==false and C:isempty()==false
        test A:dim()==2 and A:issingulardir(0)==false and A:issingulardir(1)==false
        test B:dim()==1 and B:issingulardir(0)==false and B:issingulardir(1)==true
        test C:dim()==1 and C:issingulardir(0)==true and C:issingulardir(1)==false
        test A:vol()==8 and B:vol()==2 and C:vol()==1
    end

    testset "permutation" do
        test A.perm[0] == 0 and A.perm[1] == 1
        test B.perm[0] == 0 and B.perm[1] == 1
        test C.perm[0] == 1 and C.perm[1] == 0
    end

    testset "evaluation" do
        terracode
            var x = point_2d.from(0,1)
        end
        --test from point, list, and tuple
        test A(point_2d.from(0,0)) == x and A(0,0) == x and A({0,0})==x
        test A(0.5,0.5) == point_2d.from(1,3)
        test A(1,1) == point_2d.from(2,5)
    end

    testset "square evaluation - custom origin" do
        terracode
            A:setorigin(point_2d.from(2,5))
        end
        test A(0,0) == point_2d.from(2,5)
        test A(1,1) == point_2d.from(0,1)
    end

    testset "segment evaluation - custom origin" do
        terracode
            B:setorigin(point_2d.from(2,1))
            C:setorigin(point_2d.from(0,1))
        end
        test B(0) == point_2d.from(2,1) and B(1) == point_2d.from(0,1)
        test C(0) == point_2d.from(0,1) and C(1) == point_2d.from(0,2)
    end

    testset "square - barycentric coordinate" do
        terracode
            A:setorigin(point_2d.from(2,5))
            var x_1 = A:barycentriccoord(point_2d.from(2,5))
            var x_2 = A:barycentriccoord(point_2d.from(0,1))
        end
        test x_1 == point_2d.from(0,0) and x_2 == point_2d.from(1,1)
    end

    testset "line - barycentric coordinate" do
        terracode
            var B = line.new(interval.new(0, 2), 1)
            B:setorigin(point_2d.from(2,1))
            var x_1 = B:barycentriccoord(point_2d.from(2,1))
            var x_2 = B:barycentriccoord(point_2d.from(0,1))
        end
        test x_1 == point_1d.from(0) and x_2 == point_1d.from(1)
    end

end