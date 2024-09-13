require("terralibext")
local io = terralib.includec("stdio.h")

local math = require("math")
local vec = require("vector")
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

local Interval = {}
Interval.__index = Interval;
Interval.__metatable = Interval;

Interval.new = function(a, b)

    assert(type(a)=="number" and type(b)=="number" and "a and b need to be real numbers.")
    assert(a < b)

    --new table
    local interval = {}

    --accessible values
    interval.isinterval = true
    interval.eltype = T
    interval.a = a
    interval.b = b
    interval.volume = math.abs(b-a)

    --metatable
    return setmetatable(interval, Interval)
end

function Interval:isinside(x) return (self.a <= x) and (x <= self.b) end

function Interval:__tostring()
    return "interval("..tostring(self.a)..", " ..tostring(self.b)..")"
end

function Interval:__eq(other)
    return self.a==other.a and self.b==other.b
end

function Interval.intersection(self, other)
    
    --check input arguments
    for i,v in ipairs({self, other}) do
        if not (type(v)=="number" or type(v)=="table" and v.isinterval) then
            error("Ecpected numbers and/or intervals as input.")
        end
    end

    local intersection_point_point = function(point_a, point_b)
        if point_a==point_b then
            return point_a
        else
            return nil
        end
    end

    local intersection_interval_point = function(interval, point)
        if interval:isinside(point) then
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
                return Interval.new(a, b)
            elseif a==b then
                --intersection is a point
                return a
            else
                --no common intersection
                return nil 
            end
        end
    end

    local get_intersection_type = function(a, b)
        if type(a)=="number" and type(b)=="number" then
            return intersection_point_point(a, b)
        elseif type(a)=="number" and type(b)=="table" and b.isinterval then
            return intersection_interval_point(b, a)
        elseif type(b)=="number" and type(a)=="table" and a.isinterval then
            return intersection_interval_point(a, b)
        elseif type(a)=="table" and type(b)=="table" and a.isinterval and b.isinterval then
            return intersection_interval_interval(a, b)
        else
            error("Function arguments need to be a primitives or intervals.")
        end
    end

    return get_intersection_type(self, other)
end


import "terratest/terratest"

testenv "Interval" do

    testset "static data" do
        local interval = Interval.new(0,2)
        test [interval.isinterval]
        test [interval.eltype == T]
        test [interval.a == 0]
        test [interval.b == 2]
        test [interval:isinside(0) and interval:isinside(2) and not interval:isinside(-0.01) and not interval:isinside(2.01)]
    end

    testset "equality" do
        test [Interval.new(0,2) == Interval.new(0,2)]
        test [Interval.new(0,2) ~= Interval.new(0,1)]
    end

    testset "intersection - types" do
        local I = Interval.new(0,2)
        local J = Interval.new(1,3)
        local K = Interval.new(2,4)
        local L = Interval.new(3,5)
        test [Interval.intersection(I,I)==I]
        test [Interval.intersection(I,J)==Interval.new(1,2)]
        test [Interval.intersection(I,K)==2]
        test [Interval.intersection(I,L)==nil]
    end

end


local Hypercube = {}
Hypercube.__index = Hypercube;
Hypercube.__metatable = Hypercube;

Hypercube.new = function(...)

    local args = terralib.newlist{...}
    local N = #args

    --generate inverse permutation
    local invperm = terralib.newlist()
    local volume = 1
    local D = 0
    local alpha, beta = 1, N
    for i, v in ipairs(args) do
        if type(v)=="number" then
            invperm[i] = beta
            beta = beta - 1
        elseif v.isinterval then
            invperm[i] = alpha
            D = D + 1
            volume = volume * v.volume
            alpha = alpha + 1  
        else
            error("Arguments should be an interval or a number.")
        end
    end
    --generate permutation
    local perm = terralib.newlist()
    for i = 1, N do
        perm[ invperm[i] ] = i
    end

    --new table
    local hypercube = {}

    --static data members
    hypercube.ishypercube = true
    hypercube.rangedim = N
    hypercube.dim = D
    hypercube.volume = volume
    hypercube.I = args
    hypercube.perm = perm
    hypercube.invperm = invperm

    local typename = terralib.newlist()
    for k = 1, N do
        typename:insert(tostring(hypercube.I[k]))
    end
    hypercube.typename = "hypercube(" .. table.concat(typename,",") ..")" 
    
    return setmetatable(hypercube, Hypercube)
end

function Hypercube:__tostring()
    return self.typename
end

function Hypercube:__eq(other)
    if self.rangedim==other.rangedim then
        for i=1, self.rangedim do 
            if self.I[i] ~= other.I[i] then
                return false
            end
        end
        return true
    end
    return false
end

function Hypercube:issingulardir(i)
    return self.invperm[i] > self.dim
end

Hypercube.intersection = function(...)

    --check input arguments
    local args = terralib.newlist{...}
    if #args<2 then
        error("Expected two or more hypercubes as input.")
    end
    for i,v in ipairs(args) do
        if not (type(v)=="table" and v.ishypercube) then
            error("Expected two or more hypercubes as input.")
        end
    end
    for i,v in ipairs(args) do
        if not (v.rangedim==args[1].rangedim) then
            error("Expected hypercubes with range dimension "..tostring(rangedim))
        end
    end

    --compute intersection of two hypercubes
    local get_intersection_two_vars = function(a, b)
        --if a == b then
        --    return a
        --else
            local N = a.rangedim
            local I = terralib.newlist()
            for i = 1, N do
                local intersection = Interval.intersection(a.I[i], b.I[i])
                if intersection==nil then
                    --no common intersection
                    return nil
                end
                I:insert(intersection)
            end
            return Hypercube.new(unpack(I))
        --end
    end

    --use recursion to compute intersection of multiple hypercubes
    local get_intersection

    get_intersection = function(a, b, ...)
        local args = terralib.newlist{...}
        local t = get_intersection_two_vars(a, b)
        if #args>0 then
            t = get_intersection(t, ...)
        end
        return t
    end

    return get_intersection(...)
end

Hypercube.__mul = function(self, other)
    local cube = Hypercube.intersection(self, other) --intersection type
    if cube then --non-empty intersection
        local I = cube.I
        local N = self.rangedim
        for i = 1, N do
            if self:issingulardir(i) and not other:issingulardir(i) then
                I[i] = other.I[i]
            elseif other:issingulardir(i) and not self:issingulardir(i) then
                I[i] = self.I[i]
            elseif self:issingulardir(i) and other:issingulardir(i) then
                I[i] = self.I[i]
            else
                error("This branch should be unreachable.")
            end
        end
        return Hypercube.new(unpack(I))
    end
end

Hypercube.__div = function(self, other)
    local cube = Hypercube.intersection(self, other) --intersection type
    if cube then --non-empty intersection
        local I = cube.I
        local N = self.rangedim
        for i = 1, N do
            if other:issingulardir(i) then
                I[i] = self.I[i]
            else
                I[i] = cube.I[i].a
            end
        end
        return Hypercube.new(unpack(I))
    end
end


testenv "Hypercube - 3d" do

    local I, J, K = Interval.new(0,1), Interval.new(1,3), Interval.new(2,5)

    testset "static data" do
        local Cube = Hypercube.new(I, J, K)
        test [Cube.ishypercube]
        test [Cube.dim==3]
        test [Cube.rangedim==3]
        test [Cube.volume==6]
    end

    testset "intersection" do
        --lines
        test [Hypercube.intersection(Hypercube.new(I, 0, 0), Hypercube.new(0, I, 0)) == Hypercube.new(0, 0, 0)]
        --volumes
        test [Hypercube.intersection(Hypercube.new(I, I, I), Hypercube.new(I, I, I)) == Hypercube.new(I, I, I)]
        test [Hypercube.intersection(Hypercube.new(I, I, I), Hypercube.new(J, I, I)) == Hypercube.new(1, I, I)]
        test [Hypercube.intersection(Hypercube.new(I, I, I), Hypercube.new(J, J, I)) == Hypercube.new(1, 1, I)]
        test [Hypercube.intersection(Hypercube.new(I, I, I), Hypercube.new(J, J, J)) == Hypercube.new(1, 1, 1)]
        test [Hypercube.intersection(Hypercube.new(I, I, I), Hypercube.new(K, J, J)) == nil]
    end

    testset "mul" do
        local L1, L2, L3 = Hypercube.new(I, J.a, K.a), Hypercube.new(I.a, J, K.a), Hypercube.new(I.a, J.a, K)
        test [L1*L2*L3 == Hypercube.new(I, J, K)]
    end

    testset "div" do
        local A = Hypercube.new(Interval.new(0,2), Interval.new(0,1), Interval.new(0,1))
        local B = Hypercube.new(0, Interval.new(0,1), Interval.new(0,1))
        local C = Hypercube.new(Interval.new(0,2), 0, 0)
        test [(B * C) == A]
        test [(B * C) / C == B]
        test [(B * C) / B == C]
    end

end

local Mapping = terralib.memoize(function(args)

    local domain = args.domain
    local origin = args.origin

    --check inputs
    if domain==nil or origin==nil then
        error("Expected named arguments 'domain' and origin")
    end

    local N = domain.rangedim
    
    if not (type(domain)=="table" and domain.ishypercube) then
        error("Expected named argument 'domain' to be a hypercube.")
    end
    if not (type(origin)=="table" and #origin==N) then
        error("Expected named argument 'origin' to be an array of "..N .. " numbers.")
    end

    --dummy struct
    local struct mapping{ }

    --static data
    mapping.ismapping = true
    mapping.domain = domain
    mapping.origin = origin


    --generate inverse permutation, and fill A and B with data for linear
    --map: A * X + B, such that [0,1]^N maps to the hypercube
    local A, B = terralib.newlist(), terralib.newlist()
    local A_inv, B_inv = terralib.newlist(), terralib.newlist()
    local D = 0
    for i, v in ipairs(domain.I) do
        local o = origin[i]
        if type(v)=="number" then
            if v~=o then
                error("Expected origin["..i.."] = ".. tostring(v))
            end
            A:insert(0)
            B:insert(v)
            A_inv:insert(0)
            B_inv:insert(0)
        elseif v.isinterval then
            local center = 0.5 * (v.a + v.b)
            local signed = (o<=center) and 1 or -1
            A:insert((v.b - v.a) * signed)
            B:insert(o)
            A_inv:insert(signed / (v.b - v.a))
            B_inv:insert(-signed * o / (v.b - v.a))
            D = D + 1
        end
    end

    --construct static arrays for linear map 'y = a * x + b'
    --and origin 'o'
    local map, invmap = {}, {}
    map.a = terralib.constant(terralib.new(T[N], A))
    map.b = terralib.constant(terralib.new(T[N], B))
    map.o = terralib.constant(terralib.new(T[N], origin))
    --construct static arrays for linear inverse map 'x = (1/a) y - b/a'
    invmap.a = terralib.constant(terralib.new(T[N], A_inv))
    invmap.b = terralib.constant(terralib.new(T[N], B_inv))

    local perm = domain.perm
    local invperm = domain.invperm

    mapping.metamethods.__apply = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        if #args ~= D then
            error("Expected " .. D .. " input arguments.")
        end
        return quote
            var y = [map.o]
            escape
                for i = 1, D do
                    local x = args[i]
                    local k = perm[i]
                    local a, b = A[k], B[k]
                    emit quote y[k-1] = tmath.fusedmuladd([T](a), [T](x), [T](b)) end
                end
            end
        in
            y
        end
    end))

    local point_nd = Point(T, N)

    terra mapping:barycentriccoord(y : point_nd)
        var x : ntuple(T, D)
        escape
            for i = 1, D do
                local k = perm[i]
                local sx = "_"..tostring(i-1)
                local a, b = A_inv[k], B_inv[k]
                emit quote x.[sx] = tmath.fusedmuladd([T](a), y(k-1), [T](b)) end
            end
        end
        return x
    end


    return mapping
end)

testenv "Mapping - 3d - line" do

    local Point = Point(T, 3)

    --compute hypercube at compile-time
    local Line = Hypercube.new(1, Interval.new(2,4), 2)
    
    --compute mappings at compile-time
    local F = Mapping{domain=Line, origin={1,2,2}}
    local G = Mapping{domain=Line, origin={1,4,2}}

    testset "properties" do
        test [F.ismapping]
        test [F.domain.dim==1]
        test [F.domain.rangedim==3]
    end

    terracode
        var f : F
        var g : G
    end

    testset "evaluation" do
        terracode
            var a = f(0)
            var b = f(1)
        end
        test a[0]==1 and a[1]==2 and a[2]==2
        test b[0]==1 and b[1]==4 and b[2]==2
    end

    testset "evaluation - reversed origin" do
        terracode
            var a = g(0)
            var b = g(1)
        end
        test a[0]==1 and a[1]==4 and a[2]==2
        test b[0]==1 and b[1]==2 and b[2]==2
    end

    testset "barycentric coordinates" do
        terracode
            var a = f:barycentriccoord(Point.from(1,2,2))
            var b = f:barycentriccoord(Point.from(1,4,2))
        end
        test a._0==0
        test b._0==1
    end

    testset "barycentric coordinates - reversed origin" do
        terracode
            var a = g:barycentriccoord(Point.from(1,2,2))
            var b = g:barycentriccoord(Point.from(1,4,2))
        end
        test a._0==1
        test b._0==0
    end

end

testenv "Mapping - 3d - surface" do

    local Point = Point(T, 3)

    --compute hypercube at compile-time
    local Surf = Hypercube.new(1, Interval.new(2,4), Interval.new(5,6))
    
    --compute mappings at compile-time
    local F = Mapping{domain=Surf, origin={1,2,5}}
    local G = Mapping{domain=Surf, origin={1,4,6}}

    testset "properties" do
        test [F.ismapping]
        test [F.domain.dim==2]
        test [F.domain.rangedim==3]
    end

    terracode
        var f : F
        var g : G
    end

    testset "evaluation" do
        terracode
            var a = f(0,0)
            var b = f(1,1)
        end
        test a[0]==1 and a[1]==2 and a[2]==5
        test b[0]==1 and b[1]==4 and b[2]==6
    end

    testset "evaluation - reversed origin" do
        terracode
            var a = g(0,0)
            var b = g(1,1)
        end
        test a[0]==1 and a[1]==4 and a[2]==6
        test b[0]==1 and b[1]==2 and b[2]==5
    end

    testset "barycentric coordinates" do
        terracode
            var a = f:barycentriccoord(Point.from(1,2,5))
            var b = f:barycentriccoord(Point.from(1,4,6))
        end
        test a._0==0 and a._1==0
        test b._0==1 and b._1==1
    end

    testset "barycentric coordinates - reversed origin" do
        terracode
            var a = g:barycentriccoord(Point.from(1,2,5))
            var b = g:barycentriccoord(Point.from(1,4,6))
        end
        test a._0==1 and a._1==1
        test b._0==0 and b._1==0
    end

end


testenv "Mapping - 3d - volume" do

    local Point = Point(T, 3)

    --compute hypercube at compile-time
    local Vol = Hypercube.new(Interval.new(0,1), Interval.new(2,4), Interval.new(5,6))
    
    --compute mappings at compile-time
    local F = Mapping{domain=Vol, origin={0,2,5}}
    local G = Mapping{domain=Vol, origin={1,4,6}}

    testset "properties" do
        test [F.ismapping]
        test [F.domain.dim==3]
        test [F.domain.rangedim==3]
    end

    terracode
        var f : F
        var g : G
    end

    testset "evaluation" do
        terracode
            var a = f(0,0,0)
            var b = f(1,1,1)
        end
        test a[0]==0 and a[1]==2 and a[2]==5
        test b[0]==1 and b[1]==4 and b[2]==6
    end

    testset "evaluation - reversed origin" do
        terracode
            var a = g(0,0,0)
            var b = g(1,1,1)
        end
        test a[0]==1 and a[1]==4 and a[2]==6
        test b[0]==0 and b[1]==2 and b[2]==5
    end

    testset "barycentric coordinates" do
        terracode
            var a = f:barycentriccoord(Point.from(0,2,5))
            var b = f:barycentriccoord(Point.from(1,4,6))
        end
        test a._0==0 and a._1==0 and a._2==0
        test b._0==1 and b._1==1 and b._2==1
    end

    testset "barycentric coordinates - reversed origin" do
        terracode
            var a = g:barycentriccoord(Point.from(0,2,5))
            var b = g:barycentriccoord(Point.from(1,4,6))
        end
        test a._0==1 and a._1==1 and a._2==1
        test b._0==0 and b._1==0 and b._2==0
    end

end



local Pyramid = function(args)

    local base = args.base
    local apex = args.apex

    --check inputs
    if base==nil or apex==nil then
        error("Expected named arguments 'base' and 'apex'")
    end

    local M = base.dim+1
    local N = base.rangedim
    
    if not (type(base)=="table" and base.ishypercube) then
        error("Expected named argument 'base' to be a hypercube.")
    end
    if not (type(apex)=="table" and #origin==N) then
        error("Expected named argument 'apex' to be an array of "..N .. " numbers.")
    end

    --dummy struct
    local pyramid = { }

    --static data
    pyramid.ispyramid = true
    pyramid.base = base
    mapping.apex = apex


    function pyramid:dim()
        return self.base.dim+1
    end


    function pyramid:height()
        --local x = self.base:barycentriccoords(apex)
        --local y = self.base(x)
        return 1.0
    end

    function pyramid:vol()
        --return self:height() * self.base:vol() / self:dim()
    end


end