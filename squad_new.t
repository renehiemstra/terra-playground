require("terralibext")
local io = terralib.includec("stdio.h")

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

    --local values
    local center = 0.5 * (a + b)
    local volume = math.abs(b - a)

    --accessible values / functions
    interval.isinterval = true
    interval.eltype = T
    interval.a = a
    interval.b = b
    interval.volume = volume
    interval.staticmethods = {}

    --set dispatch static and dynamic methods
    interval.metamethods.__getmethod = function(self, methodname)
        return self.methods[methodname] or interval.staticmethods[methodname]
    end

    --static 'isinside' lua function
    interval.isinside = function(x) return (a <= x) and (x <= b) end

    --dynamic 'isinside' terra function
    terra interval:isinside(x : T) return (a <= x) and (x <= b) end

    --set the origin at runtime
    terra interval:setorigin(origin : T)
        self.origin = origin
        self.s = terralib.select(origin <= center, 1., -1.)
    end

    terra interval:getorigin()
        return self.origin
    end

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

local intersection_point_point = function(a, b)
    if a==b then
        return a
    else
        return nil
    end
end

local intersection_interval_point = function(interval, point)
    if interval.isinside(point) then
        return point
    else
        return nil
    end
end

local intersection_interval_interval = function(self, other)
    if self == other then
        return self
    else
        local a, b = math.max(self.a, other.a), math.min(self.b, other.b)
        if a < b then
            --intersection is a non-empty interval
            return Interval(a, b)
        elseif a==b then
            --intersection is a point
            return a
        else
            --no common intersection
            return nil 
        end
    end
end

local intersection_interval = terralib.memoize(function(a, b)
    --try the case of 'intervals' and 'primiatives'
    if type(a)=="number" and type(b)=="number" then
        intersection_point_point(a, b)
    elseif type(a)=="number" and type(b)=="table" and b.isinterval then
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
                local intersection = intersection_interval(a.I[i], b.I[i])
                if intersection==nil then
                    --no common intersection
                    return nil
                end
                I:insert(intersection)
            end
            return Hypercube(unpack(I))
        end
    else
        error("Function arguments need to be two hypercubes.")
    end
end

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

    testset "evaluation, inside" do
        terracode
            var I : Interval(0, 2)
        end
        test I(0)==0 and I(1)==2
        terracode
            I:setorigin(2)
        end
        test I(0)==2 and I(1)==0
        test I:isinside(0.0) and I:isinside(2.0) and not I:isinside(-0.01) and not I:isinside(2.01)
    end

    testset "intersection - types" do
        local A = Interval(1,3)
        local B = Interval(2,4)
        local C = Interval(3,4)
        local D = Interval(4,5)
        test [intersection_interval(A, B) == Interval(2,3)]
        test [intersection_interval(A, D) == nil]
    end

    testset "intersection" do
        terracode
            var u : Interval(1,3)
            var v : Interval(2,4)
            var z_1 = intersect(v, u)
        end
        test z_1(0) == 2 and z_1(1) == 3
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

    --generate inverse permutation
    local invperm = terralib.newlist()
    local volume = 1
    local D = 0
    local Origin = terralib.newlist()
    local alpha, beta = 1, N
    for i, v in ipairs(args) do
        if type(v)=="number" then
            invperm[i] = beta
            beta = beta - 1
            Origin:insert(v)
        elseif v.isinterval then
            invperm[i] = alpha
            D = D + 1
            volume = volume * v.volume
            alpha = alpha + 1
            Origin:insert(0)  
        else
            error("Arguments should be an interval or a number.")
        end
    end
    --generate permutation
    local perm = terralib.newlist()
    for i = 1, N do
        perm[ invperm[i] ] = i
    end

    --definition of struct
    local hypercube = terralib.types.newstruct("hypercube")
    --add entries (only the intervals)
    for k = 1, D do
        local i = perm[k]
        local ftype = args[i]
        local fname = "_"..tostring(k-1)
        hypercube.entries:insert({field = fname, type = ftype})
    end
	hypercube:complete()

    hypercube.ishypercube = true
    hypercube.volume = volume
    hypercube.rangedim = N
    hypercube.dim = D
    hypercube.I = args

    local point_nd = Point(T, N)

    --hypercube origin in physical space
    local origin = terralib.constant(terralib.new(T[N], Origin))

    --set hypercube origin in parameter space
    hypercube.methods.setorigin = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        if #args == N then
            return quote
                escape
                    for i = 1, D do
                        local k = perm[i]
                        local y = args[k]
                        local field = "_"..tostring(i-1)
                        emit quote self.[field]:setorigin(y) end
                    end
                end
            end
        else
            error("Function requires " .. N .. " input arguments.")
        end
    end))

    function hypercube.issingulardir(i)
        return invperm[i] > D
    end

    hypercube.metamethods.__apply = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        if #args == D then
            return quote
                var y = origin
                escape
                    for i = 1, D do
                        local x = args[i]
                        local k = perm[i]
                        local field = "_"..tostring(i-1)
                        emit quote y[k-1] = self.[field](x) end
                    end
                end
            in
                y
            end
        end
    end))

    terra hypercube:barycentriccoord(y : point_nd)
        var x : ntuple(T, D)    
        escape
            for i = 1, D do
                local k = perm[i]
                local field = "_"..tostring(i-1)
                emit quote x.[field] = self.[field]:barycentriccoord(y(k-1)) end
            end
        end
        return x
    end

    return hypercube
end)



testenv "hypercube - lines" do

    local Point = Point(T, 2)
    local Line = Hypercube(2, Interval(3,4))

    terracode
        var line : Line
    end

    testset "properties" do
        test [Line.dim==1]
        test [Line.rangedim==2]
    end

    testset "evaluation - line" do
        terracode
            var a, b = line(0), line(1)
        end
        test a[0]==2 and a[1]==3
        test b[0]==2 and b[1]==4
    end

    testset "evaluation - line - reversed origin" do
        terracode
            line:setorigin(2,4)
            var a, b = line(0), line(1)
        end
        test b[0]==2 and b[1]==3
        test a[0]==2 and a[1]==4
    end

    testset "barycentric coordinates" do
        terracode
            var xa = line:barycentriccoord(Point.from(2, 3))
            var xb = line:barycentriccoord(Point.from(2, 4))
        end
        test xa._0==0 and xb._0==1 
    end

    testset "barycentric coordinates - reversed origin" do
        terracode
        line:setorigin(2,4)
            var xa = line:barycentriccoord(Point.from(2, 3))
            var xb = line:barycentriccoord(Point.from(2, 4))
        end
        test xa._0==1 and xb._0==0
    end

    testset "intersection-types" do
        local A = Hypercube(0, Interval(0,1))
        local B = Hypercube(Interval(0,1), 0)
        local C = Hypercube(Interval(0,1), 2)
        test [intersection_hypercube(A, B) == Hypercube(0, 0)]
        test [intersection_hypercube(A, C) == nil]
        test [intersection_hypercube(B, C) == nil]
    end

    testset "intersection" do
        terracode
            var l1 : Hypercube(0, Interval(0,1))
            var l2 : Hypercube(Interval(0,1), 0)
            var I = intersect(l1, l2)
            var x = I()
        end
        test x[0]==0 and x[1]==0
    end

    testset "product - types" do
        local A = Hypercube(0, Interval(0,1))
        local B = Hypercube(Interval(0,1), 0)
        test [hypercube_mult(A, B)==Hypercube(Interval(0,1), Interval(0,1))]
    end

    testset "product" do
        terracode
            var l1 : Hypercube(0, Interval(0,1))
            var l2 : Hypercube(Interval(0,1), 0)
            var s = l1 * l2 --surface
        end
        test
    end

end




testenv "hypercube - lines" do

    local point = Point(T, 2)

    terracode
        var line : Hypercube(2, Interval(3,4))
        var square : Hypercube(Interval(1,3), Interval(3,4))
        var surface : Hypercube(4, Interval(1,3), Interval(3,4))
        var volume : Hypercube(Interval(0,1), Interval(2,4), Interval(5,8))
    end

    testset "evaluation - square" do
        terracode
            var a = square(0, 0)
            var b = square(1, 1)
        end
        test a[0]==1 and a[1]==3
        test b[0]==3 and b[1]==4
    end

    testset "evaluation - surface" do
        terracode
            var a = surface(0, 0)
            var b = surface(1, 1)
            var c = surface(0, 1)
        end
        test a[0]==4 and a[1]==1 and a[2]==3
        test b[0]==4 and b[1]==3 and b[2]==4
        test c[0]==4 and c[1]==1 and c[2]==4
    end

    testset "evaluation - volume" do
        terracode
            var a = volume(0,0,0)
            var b = volume(1,1,1)
        end
        test a[0]==0 and a[1]==2 and a[2]==5
        test b[0]==1 and b[1]==4 and b[2]==8
        terracode
            volume:setorigin(1, 4, 8)
            a = volume(0,0,0)
            b = volume(1,1,1)
        end
        test b[0]==0 and b[1]==2 and b[2]==5
        test a[0]==1 and a[1]==4 and a[2]==8
    end


    testset "intersection-types - surface" do
        local A = Hypercube(0, Interval(0,1), Interval(0,1))
        local B = Hypercube(Interval(0,1), Interval(0,1), 0)
        test [intersection_hypercube(A, B) == Hypercube(0, Interval(0,1), 0)]
    end

    testset "intersection-types - volumes" do
        local A = Hypercube(Interval(0,1), Interval(0,1), Interval(0,1))
        local B = Hypercube(Interval(1,2), Interval(0,1), Interval(0,1))
        local C = Hypercube(Interval(1,2), Interval(1,2), Interval(0,1))
        local D = Hypercube(Interval(1,2), Interval(1,2), Interval(1,2))
        test [intersection_hypercube(A, B) == Hypercube(1, Interval(0,1), Interval(0,1))]
        test [intersection_hypercube(A, C) == Hypercube(1, 1, Interval(0,1))]
        test [intersection_hypercube(A, D) == Hypercube(1, 1, 1)]
    end



end