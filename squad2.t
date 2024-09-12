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

local Interval = {}


Interval.new = terralib.memoize(function(a, b)

    assert(type(a)=="number" and type(b)=="number" and "a and b need to be real numbers.")
    assert(a < b)

    local struct interval{
    }

    --accessible values / functions
    interval.isinterval = true
    interval.eltype = T
    interval.a = a
    interval.b = b

    interval.metamethods.__typename = function(self)
        return "interval("..tostring(a).."," ..tostring(b)..")"
    end

    --static 'isinside' lua function
    interval.isinside = function(x) return (a <= x) and (x <= b) end

    --dynamic 'isinside' terra function
    terra interval:isinside(x : T) return (a <= x) and (x <= b) end

    interval.metamethods.__eq = macro(function(self, other)
        local self_type, other_type = self.tree.type, other.tree.type 
        return `[self_type.a==other_type.a and self_type.b==other_type.b]
    end)

    return interval
end)

Interval.intersection_types = function(...)

    local args = terralib.newlist{...}
    if #args<2 then
        error("Expected two or more input arguments.")
    end
    for i,v in ipairs(args) do
        if not (type(v)=="number" or type(v)=="table" and v.isinterval) then
            error("Ecpected numbers and/or intervals as input.")
        end
    end

    local intersection_point_point = function(a, b)
        if a==b then
            return a
        else
            return nil
        end
    end

    local intersection_interval_point = function(self, point)
        if self.isinside(point) then
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

    local get_intersection_type_two_vars = function(a, b)
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


    local get_intersection_type

    get_intersection_type = function(a_type, b_type, ...)
        local args = terralib.newlist{...}
        local type = get_intersection_type_two_vars(a_type, b_type)
        if #args>0 then
            type = get_intersection_type(type, ...)
        end
        return type
    end

    return get_intersection_type(...)
end

Interval.intersection = macro(function(...)
    local args = terralib.newlist{...}
    if #args<2 then
        error("Expected two or more input arguments.")
    end
    local types = terralib.newlist()
    for i,v in ipairs(args) do
        types:insert(v.tree.type)
    end
    local t = Interval.intersection_types(unpack(types))
    if t then
        if type(t)=="number" then
            return `t
        else
            return quote
                var intersection_obj : t
            in
                intersection_obj
            end
        end
    end
end)

import "terratest/terratest"

testenv "interval" do

    testset "static data" do
        local interval = Interval.new(0,2)
        test [interval.isinterval]
        test [interval.eltype == T]
        test [interval.a == 0]
        test [interval.b == 2]
        test [interval.isinside(0) and interval.isinside(2) and not interval.isinside(-0.01) and not interval.isinside(2.01)]
    end

    testset "inside" do
        terracode
            var l : Interval.new(0, 2)
        end
        test l:isinside(0.0) and l:isinside(2.0) and not l:isinside(-0.01) and not l:isinside(2.01)
    end

    testset "intersection - types" do
        local I = Interval.new(0,2)
        local J = Interval.new(1,3)
        local K = Interval.new(2,4)
        local L = Interval.new(3,5)
        test [Interval.intersection_types(I,I)==I]
        test [Interval.intersection_types(I,J)==Interval.new(1,2)]
        test [Interval.intersection_types(I,K)==2]
        test [Interval.intersection_types(I,L)==nil]
        test [Interval.intersection_types(I,I,I)==I] --test multiple inputs
    end

    testset "intersection - types" do
        terracode
            var i : Interval.new(0,2)
            var j : Interval.new(1,3)
            var k : Interval.new(2,4)
            var l : Interval.new(3,5)

            var m = Interval.intersection(i,j)
        end
        test Interval.intersection(i,i)==i
        test Interval.intersection(i,j)==m
        test Interval.intersection(i,k)==2
        test Interval.intersection(i,i,i)==i --test multiple inputs
    end

end


