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

    interval.metamethods.__typename = function(self)
        return "interval("..tostring(a).."," ..tostring(b)..")"
    end

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

local intersection_interval = function(a, b)
    --try the case of 'intervals' and 'primiatives'
    if type(a)=="number" and type(b)=="number" then
        return intersection_point_point(a, b)
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
end

local Hypercube 

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

local get_intersection_type_two_vars = function(a_type, b_type)
    local type = nil
    if a_type.ishypercube and b_type.ishypercube then
        --try hypercubes
        type = intersection_hypercube(a_type, b_type)
    else
        --try intervals
        type = intersection_interval(a_type, b_type)
    end
    return type
end

local get_intersection_type

get_intersection_type = function(a_type, b_type, ...)
    local args = terralib.newlist{...}
    local type = get_intersection_type_two_vars(a_type, b_type)
    if #args>0 then
        type = get_intersection_type(type, ...)
    end
    return type
end

local intersect = macro(function(...)
    local args = terralib.newlist{...}
    if #args<2 then
        error("Expected two or more input arguments.")
    end
    local types = terralib.newlist()
    for i,v in ipairs(args) do
        types:insert(v.tree.type)
    end
    local type = get_intersection_type(unpack(types))
    if type then
        return quote
            var intersection_obj : type
        in
            intersection_obj
        end
    end
end)

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
    hypercube.ishypercube = true
    hypercube.volume = volume
    hypercube.rangedim = N
    hypercube.dim = D
    hypercube.I = args

    --add entries
    local typename = terralib.newlist()
    for k = 1, N do
        local fname, ftype = "_"..tostring(k-1), nil
        if type(hypercube.I[k])=="number" then
            ftype = T
        elseif hypercube.I[k].isinterval then
            ftype = hypercube.I[k]
        else
            error("Expected an interval or a real number.")
        end
        hypercube.entries:insert({field = fname, type = ftype})
        typename:insert(tostring(hypercube.I[k]))
    end
	hypercube:complete()

    hypercube.metamethods.__typename = function(self)
        return "hypercube(" .. table.concat(typename,",") ..")"
    end

    --hypercube origin in physical space
    local origin = terralib.constant(terralib.new(T[N], Origin))

    --set singular dimensions
    hypercube.methods.__init = terra(self : &hypercube)
        escape
            for k = 1, N do
                local s = "_"..tostring(k-1)
                if type(hypercube.I[k])=="number" then
                    emit quote self.[s] = origin[k-1] end
                else
                    emit quote self.[s]:__init() end
                end
            end
        end
    end

    local point_nd = Point(T, N)

    --set hypercube origin in parameter space
    hypercube.methods.setorigin = macro(function(self, ...)
        local args = terralib.newlist{...}
        if #args ~= N then
            error("Expected " .. N .. " input arguments.")
        end
        return quote
            escape
                local tp = self.tree.type
                for k,v in ipairs(tp:getentries()) do
                    if v and v.type.isinterval then
                        local s = v.field
                        local y = args[k]
                        emit quote self.[s]:setorigin([y]) end
                    end
                end
            end
        end
    end)

    function hypercube.issingulardir(i)
        return invperm[i] > D
    end

    hypercube.metamethods.__apply = macro(terralib.memoize(function(self, ...)
        local args = terralib.newlist{...}
        if #args ~= D then
            error("Expected " .. D .. " input arguments.")
        end
        return quote
            var y = origin
            escape
                for i = 1, D do
                    local x = args[i]
                    local k = perm[i]
                    local s = "_"..tostring(k-1)
                    emit quote y[k-1] = self.[s](x) end
                end
            end
        in
            y
        end
    end))

    terra hypercube:barycentriccoord(y : point_nd)
        var x : ntuple(T, D)    
        escape
            for i = 1, D do
                local k = perm[i]
                local sx = "_"..tostring(i-1)
                local sy = "_"..tostring(k-1)
                emit quote x.[sx] = self.[sy]:barycentriccoord(y(k-1)) end
            end
        end
        return x
    end

    hypercube.metamethods.__mul = macro(function(self, other)
        local stype, otype = self.tree.type, other.tree.type
        local multtype = hypercube_mult(stype, otype)
        if multtype then
            return quote
                var res : multtype
                escape
                    for k = 1, N do
                        local s = "_"..tostring(k-1)
                        if stype.issingulardir(k) and not otype.issingulardir(k) then
                            emit quote res.[s]:setorigin(other.[s]:getorigin()) end
                        elseif otype.issingulardir(k) and not stype.issingulardir(k) then
                            emit quote res.[s]:setorigin(self.[s]:getorigin()) end
                        end
                    end
                end
            in
                res
            end
        end
    end)

    hypercube.metamethods.__div = macro(function(self, other)
        local stype, otype = self.tree.type, other.tree.type
        local divtype = hypercube_div(stype, otype)
        if divtype then
            return quote
                var res : divtype
                escape
                    for k,v in ipairs(divtype:getentries()) do
                        local s = v.field
                        if otype.issingulardir(k) and not stype.issingulardir(k) then
                            emit quote res.[s]:setorigin(other.[s]) end
                        end
                    end
                end
            in
                res
            end
        end
    end)

    return hypercube
end)



testenv "hypercube - 2d" do

    local Point = Point(T, 2)
    local Line = Hypercube(2, Interval(3,4))

    terracode
        var line : Line
    end

    testset "properties" do
        test [Line.ishypercube]
        test [Line.dim==1]
        test [Line.rangedim==2]
        test [Line.volume==1]
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
        local B = Hypercube(Interval(0,2), 0)
        test [hypercube_mult(A, B)==Hypercube(Interval(0,2), Interval(0,1))]
        test [hypercube_mult(A, B)==hypercube_mult(B, A)]
    end

    testset "product" do
        terracode
            var l1 : Hypercube(Interval(0,2), 0)
            var l2 : Hypercube(0, Interval(0,1))
            var surf = l1 * l2 --surface
            var a, b = surf(0,0), surf(1,1)
        end
        test a[0]==0 and a[1]==0
        test b[0]==2 and b[1]==1
    end

    testset "product - reversed origin" do
        terracode
            var l1 : Hypercube(Interval(0,2), 0)
            l1:setorigin(2, 0)
            var l2 : Hypercube(0, Interval(0,1))
            l2:setorigin(0, 1)
            var surf = l1 * l2 --surface
            var a, b = surf(0,0), surf(1,1)
        end
        test a[0]==2 and a[1]==1
        test b[0]==0 and b[1]==0
    end

    testset "div - types" do
        local A = Hypercube(Interval(0,2), Interval(0,1))
        local B = Hypercube(Interval(0,2), 0)
        local C = Hypercube(0, Interval(0,1))
        test [hypercube_div(A, B)==C]
        test [hypercube_div(A, C)==B]
    end

    testset "div" do
        terracode
            var surf : Hypercube(Interval(0,2), Interval(0,1))
            var l1 : Hypercube(Interval(0,2), 0)
            var l2 : Hypercube(0, Interval(0,1))
            var g = surf / l1
            var g_0, g_1 = g(0), g(1)
            var h = surf / l2
            var h_0, h_1 = h(0), h(1)
        end
        test g_0[0]==0 and g_0[1]==0 and g_1[0]==0 and g_1[1]==1
        test h_0[0]==0 and h_0[1]==0 and h_1[0]==2 and h_1[1]==0
    end

    testset "div - reversed origin" do
        terracode
            var surf : Hypercube(Interval(0,2), Interval(0,1))
            var l1 : Hypercube(Interval(0,2), 1)
            var l2 : Hypercube(2, Interval(0,1))
            var g = surf / l1
            var g_0, g_1 = g(0), g(1)
            var h = surf / l2
            var h_0, h_1 = h(0), h(1)
        end
        test g_0[0]==0 and g_0[1]==1 and g_1[0]==0 and g_1[1]==0
        test h_0[0]==2 and h_0[1]==0 and h_1[0]==0 and h_1[1]==0
    end

end

testenv "hypercube - 3d" do

    local Point = Point(T, 3)
    local Vol = Hypercube(Interval(0,1), Interval(2,4), Interval(5,6))
    local I, J, K = Interval(0,1), Interval(1,2), Interval(2,3)

    terracode
        var vol : Vol
    end

    testset "properties" do
        test [Vol.ishypercube]
        test [Vol.dim==3]
        test [Vol.rangedim==3]
        test [Vol.volume==2]
    end

    testset "evaluation - volume" do
        terracode
            var a = vol(0,0,0)
            var b = vol(1,1,1)
        end
        test a[0]==0 and a[1]==2 and a[2]==5
        test b[0]==1 and b[1]==4 and b[2]==6
    end

    testset "evaluation - volume - reversed origin" do
        terracode
            vol:setorigin(1, 4, 6)
            var a = vol(0,0,0)
            var b = vol(1,1,1)
        end
        test a[0]==1 and a[1]==4 and a[2]==6
        test b[0]==0 and b[1]==2 and b[2]==5
    end

    testset "barycentric coordinates" do
        terracode
            var a = vol:barycentriccoord(Point.from(0,2,5))
            var b = vol:barycentriccoord(Point.from(1,4,6))
        end
        test a._0==0 and a._1==0 and a._2==0
        test b._0==1 and b._1==1 and b._2==1
    end

    testset "barycentric coordinates - reversed origin" do
        terracode
            vol:setorigin(1, 4, 6)
            var a = vol:barycentriccoord(Point.from(0,2,5))
            var b = vol:barycentriccoord(Point.from(1,4,6))
        end
        test a._0==1 and a._1==1 and a._2==1
        test b._0==0 and b._1==0 and b._2==0
    end

    testset "intersection-types" do
        --lines
        test [intersection_hypercube(Hypercube(I, 0, 0), Hypercube(0, I, 0)) == Hypercube(0, 0, 0)]
        --volumes
        test [intersection_hypercube(Hypercube(I, I, I), Hypercube(I, I, I)) == Hypercube(I, I, I)]
        test [intersection_hypercube(Hypercube(I, I, I), Hypercube(J, I, I)) == Hypercube(1, I, I)]
        test [intersection_hypercube(Hypercube(I, I, I), Hypercube(J, J, I)) == Hypercube(1, 1, I)]
        test [intersection_hypercube(Hypercube(I, I, I), Hypercube(J, J, J)) == Hypercube(1, 1, 1)]
        test [intersection_hypercube(Hypercube(I, I, I), Hypercube(K, J, J)) == nil]
    end

    testset "intersection - 3 lines" do
        terracode
            var l1 : Hypercube(I, 0, 0)
            var l2 : Hypercube(0, I, 0)
            var l3 : Hypercube(0, 0, I)
            var p = intersect(l1, l2, l3)
            var x = p()
        end
        test x[0]==0 and x[1]==0 and x[2]==0
    end

    testset "intersection - 2 volumes" do
        terracode
            var v1 : Hypercube(I, I, I)
            var v2 : Hypercube(J, I, I)
            var s = intersect(v1, v2)
            var a, b = s(0, 0), s(1, 1)
        end
        test a[0]==1 and a[1]==0 and a[2]==0
        test b[0]==1 and b[1]==1 and b[2]==1
    end

    testset "intersection - 3 volumes" do
        terracode
            var v1 : Hypercube(I, I, I)
            var v2 : Hypercube(J, I, I)
            var v3 : Hypercube(J, J, I)
            var l = intersect(v1, v2, v3)
            var a, b = l(0), l(1)
        end
        test a[0]==1 and a[1]==1 and a[2]==0
        test b[0]==1 and b[1]==1 and b[2]==1
    end

    testset "product - types" do
        local L1, L2, L3 = Hypercube(I, 0, 0), Hypercube(0, I, 0), Hypercube(0, 0, I)
        test [hypercube_mult(L1, hypercube_mult(L3, L2))==hypercube_mult(L3, hypercube_mult(L2, L1))]
    end

    testset "product" do
        terracode
            var l1 : Hypercube(Interval(0,1), 1, 2)
            var l2 : Hypercube(0, Interval(1,2), 2)
            var l3 : Hypercube(0, 1, Interval(2,3))

            var vol = l1 * l2 * l3 --volume
            var a, b = vol(0,0,0), vol(1,1,1)
        end
        test a[0]==0 and a[1]==1 and a[2]==2
        test b[0]==1 and b[1]==2 and b[2]==3
    end

    testset "product - reversed origin" do
        terracode
            var l1 : Hypercube(Interval(0,1), 1, 2)
            l1:setorigin(1, 1, 2)
            var l2 : Hypercube(0, Interval(1,2), 2)
            l2:setorigin(0, 2, 2)
            var l3 : Hypercube(0, 1, Interval(2,3))
            l3:setorigin(0, 1, 3)

            var vol = l1 * l2 * l3 --volume
            var a, b = vol(0,0,0), vol(1,1,1)
        end
        test b[0]==0 and b[1]==1 and b[2]==2
        test a[0]==1 and a[1]==2 and a[2]==3
    end

    testset "div - types" do
        local A = Hypercube(Interval(0,2), Interval(0,1), Interval(0,1))
        local B = Hypercube(0, Interval(0,1), Interval(0,1))
        local C = Hypercube(Interval(0,2), 0, 0)
        local X1, X2 = hypercube_div(A, B), hypercube_div(A, C)
        test [X1==C]
        test [X2==B]
    end

    testset "div" do
        terracode
            var vol : Hypercube(Interval(0,2), Interval(0,1), Interval(0,1))
            var surf : Hypercube(2, Interval(0,1), Interval(0,1))
            var line : Hypercube(Interval(0,2), 1, 1)
            var g = vol / surf
            var g0, g1 = g(0), g(1)
            var h = vol / line
            var h0, h1 = h(0,0), h(1,1)
        end
        test g0[0]==2 and g0[1]==0 and g0[2]==0 
        test g1[0]==0 and g1[1]==0 and g1[2]==0
        test h0[0]==0 and h0[1]==1 and h0[2]==1 
        test h1[0]==0 and h1[1]==0 and h1[2]==0
        
    end
end
