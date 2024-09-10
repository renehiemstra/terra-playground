require("terralibext")

local io = terralib.includec("stdio.h")
local fun = require("fun")

local math = require("math")
local tmath = require("mathfuns")
local svector = require("svector")
local err = require("assert")

local size_t = uint64
local T = double

--return a terra tuple type of length N: {T, T, ..., T}
local ntuple = function(T, N)
    local t = terralib.newlist()
    for i = 1, N do
        t:insert(T)
    end
    return tuple(unpack(t))
end

local Point = svector.StaticVector
local Hypercube

local Interval = terralib.memoize(function(a, b)

    assert(type(a)=="number" and type(b)=="number" and "a and b need to be real numbers.")
    assert(a < b)

    local struct interval{
        origin : T
        s : int8
    }

    terra interval:__init()
        self.origin = a
        self.s = 1
    end

    interval.staticmethods = {}

    --local values
    local center = 0.5 * (a + b)
    local volume = math.abs(b - a)

    --accessible values / functions
    interval.isinterval = true
    interval.eltype = T
    interval.a = a
    interval.b = b
    interval.volume = volume

    interval.isinside = function(x) return (a <= x) and (x <= b) end

    terra interval:isinside(x : T) return (a <= x) and (x <= b) end

    interval.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or interval.staticmethods[methodname]
    end

    terra interval:setorigin(origin : T)
        self.origin = origin
        self.s = terralib.select(origin <= center, 1., -1.)
    end

    terra interval:getorigin()
        return self.origin
    end

    interval.staticmethods.new = macro(function(origin)
        assert(origin.tree.type:isprimitive())
        return quote
            var I : interval
            I:setorigin(origin)
        in
            I
        end
    end)

    --parameterization: y = A * x + o
    interval.metamethods.__apply = terra(self : interval, x : T) : T
        return tmath.fusedmuladd(x, [b-a]*self.s, self.origin)
    end

    --inverse parameterization: x = (1/A) * y - o / A
    terra interval:barycentriccoord(y : T) : T
        return tmath.fusedmuladd(y, [1./(b-a)]*self.s, [-1./(b-a)]*self.origin*self.s)
    end

    --intervals are statically determined. So objects of the same type are equal.
    interval.metamethods.__eq = terra(self : interval, other : interval) : bool
        return true
    end

    --intervals are statically determined. So objects of the same type are equal.
    interval.metamethods.__ne = terra(self : interval, other : interval) : bool
        return false
    end
    
    terra interval:vol()
        return volume
    end

    return interval
end)

local intersection

local intersection_interval_point = function(interval, point)
    if interval.isinside(point) then
        return point
    end
end

local intersection_interval_interval = function(self, other)
    if self == other then
        return self
    else
        local a, b = math.max(self.a, other.a), math.min(self.b, other.b)
        if a < b then
            --intersection is a non-empty interval
            return Interval(a, b, a)
        elseif a==b then
            --intersection is a point
            return a
        end
    end
end

local intersection_interval = terralib.memoize(function(a, b)
    --try the case of 'intervals' and 'primiatives'
    if type(a)=="number" and type(b)=="table" and b.isinterval then
        return intersection_interval_point(b, a)
    elseif type(b)=="number" and type(a)=="table" and a.isinterval then
        return intersection_interval_point(a, b)
    --try the case of two intervals
    elseif type(a)=="table" and type(b)=="table" and a.isinterval and b.isinterval then
        return intersection_interval_interval(a, b)
    else
        error("Function arguments need to be a primitives or intervals.")
    end
end)

local intersection_hypercube = function(a, b)
    if a.ishypercube and b.ishypercube then
        assert(a.rangedim == b.rangedim and "hypercubes need to have equal range dimension.")
        if a == b then
            return a
        else
            local N = a.rangedim
            local I = terralib.newlist()
            for i = 1, N do
                I:insert(intersection_interval(a.I[i], b.I[i]))
            end
            return Hypercube(unpack(I))
        end
    else
        error("Function arguments need to be two hypercubes.")
    end
end

local hypercube_mult = function(a, b)
    local cube = intersection_hypercube(a, b) --intersection type
    if cube then --non-empty intersection
        local I = cube.I
        local N = a.rangedim
        for i = 1, N do
            if a.issingulardir(i) and not b.issingulardir(i) then
                I[i] = b.I[i]
            elseif b.issingulardir(i) and not a.issingulardir(i) then
                I[i] = a.I[i]
            elseif a.issingulardir(i) and b.issingulardir(i) then
                I[i] = a.I[i]
            else
                error("This branch should be unreachable.")
            end
        end
        return Hypercube(unpack(I))
    end
end

local hypercube_div = function(a, b)
    local cube = intersection_hypercube(a, b) --intersection type
    if cube then --non-empty intersection
        local I = cube.I
        local N = a.rangedim
        for i = 1, N do
            if b.issingulardir(i) then
                I[i] = a.I[i]
            else
                I[i] = cube.I[i].a
            end
        end
        return Hypercube(unpack(I))
    end
end

Hypercube = terralib.memoize(function(...)

    local args = terralib.newlist{...}
    local N = #args
    local A, B, O = terralib.newlist(), terralib.newlist(), terralib.newlist()
    local A_inv, B_A_inv = terralib.newlist(), terralib.newlist()

    --generate inverse permutation, and fill A and B with data for linear
    --map: A * X + B, such that [0,1]^N maps to the hypercube
    local invperm = terralib.newlist()
    local volume = 1
    local D = 0
    local alpha, beta = 1, N
    for i, v in ipairs(args) do
        if type(v)=="number" then
            invperm[i] = beta
            beta = beta - 1
            A:insert(0)
            B:insert(v)
            A_inv:insert(0)
            B_A_inv:insert(0)
        elseif v.isinterval then
            invperm[i] = alpha
            A:insert(v.b - v.a)
            B:insert(v.a)
            A_inv:insert(1.0 / (v.b - v.a))
            B_A_inv:insert(v.a / (v.b - v.a))
            D = D + 1
            volume = volume * v.volume
            alpha = alpha + 1
        else
            error("Arguments should be an interval or a number.")
        end
        --fill the zero vector
        O:insert(0)
    end
    --generate permutation
    local perm = terralib.newlist()
    for i = 1, N do
        perm[ invperm[i] ] = i
    end

    local svecn = vector(T,N)
    local svecd = vector(T,D)
    local point_nd = Point(T,N)
    local point_dd = Point(T,D)
    local dtuple = ntuple(T, D)

    --construct static arrays for linear map 'y = a * x + b'
    --and origin 'o'
    local a = terralib.constant(terralib.new(T[N], A))
    local b = terralib.constant(terralib.new(T[N], B))
    local o = terralib.constant(terralib.new(T[N], O))
    --construct static arrays for linear inverse map 'x = (1/a) y - b/a'
    local a_inv = terralib.constant(terralib.new(T[N], A_inv))
    local b_a_inv = terralib.constant(terralib.new(T[N], B_A_inv))

    --dummy struct we attach the methods and static data to
    local struct hypercube{
    }
    hypercube.ishypercube = true
    hypercube.volume = volume
    hypercube.rangedim = N
    hypercube.dim = D
    hypercube.I = args

    function hypercube.issingulardir(i)
        return invperm[i] > D
    end

    local terra eval(self : hypercube, a : &svecn, b : &svecn, x : &svecn)
        return @a * @x + @b
    end

    hypercube.metamethods.__apply = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        --first try a point or a tuple
        if #args==1 and args[1].tree.type==point_dd then
            local v = args[1]
            return quote
                var a, b, x = [&svecn](&a), [&svecn](&b), o
                escape
                    for i = 1, D do
                        local k = perm[i]
                        emit quote x[k-1] = v(i-1) end
                    end
                end
            in
                eval(self, a, b, [&svecn](&x))
            end
        else
            assert(#args == D and "Number of arguments inconsistent with dimensions.")
            return quote
                var a, b, x = [&svecn](&a), [&svecn](&b), o
                escape
                    for i = 1, D do
                        local k = perm[i]
                        emit quote x[k-1] = [ args[i] ] end
                    end
                end
            in
                eval(self, a, b, [&svecn](&x))
            end
        end
    end))

    local terra barycentriccoord(a_inv : &svecn, b_a_inv : &svecn, y : &svecn)
        var x = @a_inv * @y - @b_a_inv
        var t : dtuple
        escape
            for i = 1, D do
                local k = perm[i]
                emit quote t.["_"..tostring(i-1)] = x[k-1] end
            end
        end
        return t
    end

    terra hypercube:barycentriccoord(y : svecn)
        return barycentriccoord([&svecn](&a_inv), [&svecn](&b_a_inv), &y)
    end

    terra hypercube:getorigin()
        var o : point_nd
        escape
            for i = 1, N do
                local I = hypercube.I[i]
                if type(I)=="number" then
                    emit quote o(i-1) = I end
                elseif type(I)=="table" and I.isinterval then
                    emit quote o(i-1) = I.origin end
                else
                    error("Clause should be unreachable.")
                end 
            end
        end
    end


    hypercube.metamethods.__mul = macro(function(self, other)
        local ta, tb = self.tree.type, other.tree.type
        if ta.ishypercube and tb.ishypercube then
            local type = hypercube_mult(ta, tb)
            return quote
                var I : type
            in
                I
            end
        end
    end)

    hypercube.metamethods.__div = macro(function(self, other)
        local ta, tb = self.tree.type, other.tree.type
        if ta.ishypercube and tb.ishypercube then
            local type = hypercube_div(ta, tb)
            return quote
                var I : type
            in
                I
            end
        end
    end)

    return hypercube
end)


local intersect = macro(terralib.memoize(function(self, other)
    local ta, tb = self.tree.type, other.tree.type
    local type = nil
    if ta.ishypercube and tb.ishypercube then
        --try hypercubes
        type = intersection_hypercube(ta, tb)
    else
        --try intervals
        type = intersection_interval(ta, tb)
    end
    if type then
        return quote
            var I : type
        in
            I
        end
    end
end))

local io = terralib.includec("stdio.h")

local tmath = require("mathfuns")
local err = require("assert")

import "terratest/terratest"


testenv "interval" do

    testset "static data" do
        local interval = Interval(0,2)
        test [interval.isinterval]
        test [interval.eltype == T]
        test [interval.volume == 2]
        test [interval.a == 0]
        test [interval.b == 2]
        test [interval.isinside(0) and interval.isinside(2) and not interval.isinside(-0.01) and not interval.isinside(2.01)]
    end

    testset "inside, apply" do
        terracode
            var I : Interval(0, 2)
            I:setorigin(2)
            var J = [Interval(0, 2)].new(2)
        end
        test I(0)==2 and I(1)==0
        test J(0)==2 and J(1)==0
    end

    testset "intersection" do
        terracode
            var u : Interval(1,3)
            var v : Interval(2,4)
            var w = intersect(v, u)
        end
        test w(0) == 2 and w(1) == 3
    end

    testset "barycentric coordinate" do
        terracode
            var u : Interval(1,3)
            var v : Interval(1,3)
            v:setorigin(3)
        end
        test u:barycentriccoord(1)==0 and u:barycentriccoord(3)==1
        test v:barycentriccoord(1)==1 and v:barycentriccoord(3)==0
    end

end


testenv "hypercube - lines" do

    local point = Point(T, 2)

    testset "line" do
        terracode
            var l : Hypercube(2, Interval(3,4))
            var y = l(1)
        end
        test y[0]==2 and y[1]==4
    end

    testset "square" do
        terracode
            var s : Hypercube(Interval(1,3), Interval(3,4))
            var y_00 = s(0, 0)
            var y_11 = s(1, 1)
            var p = s(point.from(1,1))
        end
        test y_00[0]==1 and y_00[1]==3
        test y_11[0]==3 and y_11[1]==4
        test p[0]==3 and p[1]==4
    end

    testset "surface" do
        terracode
            var g : Hypercube(4, Interval(1,3), Interval(3,4))
            var y_10 = g(1, 0)
            var p_10 = g(point.from(1,0))
            var x = g:barycentriccoord(vectorof(T,4,1,4))
        end
        test y_10[0]==4 and y_10[1]==3 and y_10[2]==3
        test p_10[0]==4 and p_10[1]==3 and p_10[2]==3
        test x._0==0 and x._1==1
    end

    testset "volume" do
        terracode
            var A : Hypercube(Interval(0,1), Interval(0,1), Interval(0,1))
            var B : Hypercube(Interval(1,2), Interval(0,1), Interval(0,1))
            var h = intersect(A, B)
            var y_10 = h(1, 0)
            var x = h:barycentriccoord(vectorof(T,1,0,0))
        end
        test y_10[0]==1 and y_10[1]==1 and y_10[2]==0
        test x._0==0 and x._1==0
    end

    testset "product types - 1" do
        local A = Hypercube(Interval(0,1), 1, 1)
        local B = Hypercube(0, Interval(0,1), Interval(0,1))
        local P = hypercube_mult(A, B)
        test [P == Hypercube(Interval(0,1), Interval(0,1), Interval(0,1))]
    end

    testset "product types - 2" do
        local A = Hypercube(Interval(0,1), 1, 1)
        local B = Hypercube(1, Interval(0,1), Interval(0,1))
        local P = hypercube_mult(A, B)
        test [P == Hypercube(Interval(0,1), Interval(0,1), Interval(0,1))]
    end
end



--    testset "product" do
--        terracode
--            var X : Hypercube(Interval(0,1), 1, 1)
--            var Y : Hypercube(1, Interval(0,1), Interval(0,1))
--            var f = X * Y
--            var y_111 = f(1,1,1)
--        end
--        test y_111[0]==0 and y_111[1]==0 and y_111[2]==0
--    end
--
--    testset "division" do
--        terracode
--            var Y : Hypercube(1, Interval(0,1), Interval(0,1))
--            var Z : Hypercube(Interval(0,1), Interval(0,1), Interval(0,1))
--            var r = Z / Y
--            var y_0, y_1 = r(0), r(1)
--        end
--        test y_0[0]==1 and y_0[1]==0 and y_0[2]==0
--       test y_1[0]==0 and y_1[1]==0 and y_1[2]==0
--    end