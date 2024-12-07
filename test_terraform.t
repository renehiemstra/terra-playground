-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

import "terraform"

local io = terralib.includec("stdio.h")
local base = require("base")
local concept = require("concept")

import "terratest/terratest"

local Real = concept.Real
local Integer = concept.Integer
local Float = concept.Float
local size_t = uint64


testenv "terraforming free functions" do

    testset "concrete types" do
        terraform foo(a : double, b : double)
            return a * b + 2
        end
        terraform foo(a : float, b : double)
            return a * b + 3
        end
        test foo(1., 2.) == 4.0
        test foo(float(1), 2.) == 5.0
    end
    
    testset "global method dispatch parametric concrete types" do
        terraform foo(a : T, b : T) where {T : Real}
            return a * b
        end
        terraform foo(a : int, b : T) where {T : Real}
            return a * b + 1
        end
        terraform foo(a : T, b : int) where {T : Real}
            return a * b + 2
        end
        terraform foo(a : int, b : int)
            return a * b + 3
        end
        --typical use-case in dispatch
        test foo(1., 2.) == 2.0
        test foo(1, 2.) == 3.0
        test foo(1., 2) == 4.0
        test foo(1, 2) == 5.0
    end

    testset "lobal method dispatch parametric concrete types" do
        local terraform foo(a : T, b : T) where {T : Real}
            return a * b
        end
        local terraform foo(a : int, b : T) where {T : Real}
            return a * b + 1
        end
        local terraform foo(a : T, b : int) where {T : Real}
            return a * b + 2
        end
        local terraform foo(a : int, b : int)
            return a * b + 3
        end
        --typical use-case in dispatch
        test foo(1., 2.) == 2.0
        test foo(1, 2.) == 3.0
        test foo(1., 2) == 4.0
        test foo(1, 2) == 5.0
    end

    testset "global method dispatch parametric reference types" do
        terraform foo(a : &T, b : &T) where {T : Real}
            return @a * @b
        end
        terraform foo(a : &int, b : &T) where {T : Real}
            return @a * @b + 1
        end
        terraform foo(a : &T, b : &int) where {T : Real}
            return @a * @b + 2
        end
        terraform foo(a : &int, b : &int)
            return @a * @b + 3
        end
        terracode
            var x = 1.
            var y = 2.
            var i = 1
            var j = 2
        end
        test foo(&x, &y) == 2.0
        test foo(&i, &y) == 3.0
        test foo(&x, &j) == 4.0
        test foo(&i, &j) == 5.0
    end

    testset "global method dispatch mixed parametric concrete/reference types" do
        terraform foo(a : T, b : &T) where {T : Real}
            return a + @b
        end
        terraform foo(a : &T, b : T) where {T : Real}
            return @a + b + 1
        end
        terracode
            var x = 1.
            var y = 2.
        end
        test foo(x, &y) == 3.0
        test foo(&x, y) == 4.0
    end

    testset "local method dispatch parametric reference types" do
        local terraform foo(a : &T, b : &T) where {T : Real}
            return @a * @b
        end
        local terraform foo(a : &int, b : &T) where {T : Real}
            return @a * @b + 1
        end
        local terraform foo(a : &T, b : &int) where {T : Real}
            return @a * @b + 2
        end
        local terraform foo(a : &int, b : &int)
            return @a * @b + 3
        end
        terracode
            var x = 1.
            var y = 2.
            var i = 1
            var j = 2
        end
        test foo(&x, &y) == 2.0
        test foo(&i, &y) == 3.0
        test foo(&x, &j) == 4.0
        test foo(&i, &j) == 5.0
    end

    testset "nested reference types" do
        terraform foo(a : &&T, b : &&T) where {T : Real}
            return @(@a) * @(@b)
        end
        terraform foo(a : &&int, b : &&T) where {T : Real}
            return @(@a) * @(@b) + 1
        end
        terraform foo(a : &&T, b : &&int) where {T : Real}
            return @(@a) * @(@b) + 2
        end
        terraform foo(a : &&int, b : &&int)
            return @(@a) * @(@b) + 3
        end
        terracode
            var x, y, i, j = 1., 2., 1, 2
            var rx, ry, ri, rj = &x, &y, &i, &j
        end
        test foo(&rx, &ry) == 2.0
        test foo(&ri, &ry) == 3.0
        test foo(&rx, &rj) == 4.0
        test foo(&ri, &rj) == 5.0
    end

    testset "duck-typing" do
        terraform foo(a, b)
            return a * b + 1
        end
        terraform foo(a, b, c : T) where {T : Integer}
            return a * b + c
        end
        test foo(2, 3) == 7
        test foo(2, 3, 2) == 8
    end

    testset "varargs" do
        terraform foo(a : int)
            return a + 1
        end
        terraform foo(args ...)
            var res = 1
            escape
                for i,v in ipairs(args.type.entries) do
                    local s = "_" .. tostring(i-1)
                    emit quote
                        res = res * args.[s] 
                    end
                end
            end
            return res
        end
        test foo(2) == 3
        test foo(2, 3) == 6
        test foo(2, 3, 4) == 24
    end

    testset "nearly ambiguous calls" do
        --foo<T>
        terraform foo(a : T, b : T) where {T : Float}
            return a * b
        end
        --foo<T1,T2>
        terraform foo(a : T1, b : T2) where {T1 : Float, T2 : Float}
            return a * b + 1
        end
        --both of the following calls satisfy the concept dispatch 
        --checks on all arguments. 
        test foo(1., 2.) == 2.0 --calling foo<T> and foo<T1,T2> are both valid. However,
        --foo<T> is considered to be more specialized, so we pick this method.
        test foo(float(1), 2.) == 3.0 --calling foo<T> would lead to a cast, which is 
        --not allowed, so we call foo<T1,T2>
    end

    testset "reference to local typealiases" do
        terraform foo(a : size_t)
            return a + 2
        end
        test foo(size_t(2)) == 4.0
    end    

    testset "nested namespaces" do
        local ns = {}
        ns.Float = Float
        ns.bar = {}
        ns.bar.Float = Float
        --foo<T>
        terraform ns.bar.foo(a : T, b : T) where {T : ns.bar.Float}
            return a * b
        end
        --foo<T1,T2>
        terraform ns.bar.foo(a : T1, b : T2) where {T1 : ns.Float, T2 : ns.bar.Float}
            return a * b + 1
        end
        test ns.bar.foo(1., 2.) == 2.0
        test ns.bar.foo(float(1), 2.) == 3.0 
    end

    testset "Access to parametric types in escape" do
        terraform foo(a : T1, b : T2) where {T1 : Float, T2 : Float}
            escape
                assert(T1 == double and T2 == float)
            end
            return a * b + 1
        end
        test foo(2.0, float(3.0)) == 7
    end

end


testenv "terraforming class methods" do

    --dummy struct
    local struct bar{
        index : int
    }

    terracode
        var mybar = bar{1}
    end

    testset "duck-typing" do
        terraform bar:foo(a, b)
            return a * b + self.index
        end
        terraform bar:foo(a, b, c : T) where {T : Integer}
            return a * b + c + self.index
        end
        test mybar:foo(2, 3) == 7
        test mybar:foo(2, 3, 2) == 9
    end

    testset "varargs" do
        terraform bar:foo(a : int)
            return a + self.index
        end
        terraform bar:foo(args ...)
            var res = self.index
            escape
                for i,v in ipairs(args.type.entries) do
                    local s = "_" .. tostring(i-1)
                    emit quote
                        res = res * args.[s] 
                    end
                end
            end
            return res
        end
        test mybar:foo(2) == 3
        test mybar:foo(2, 3) == 6
        test mybar:foo(2, 3, 4) == 24
    end

    testset "static methods" do
        terraform bar.sin(a : T) where {T : Real}
            return a + 2
        end
        terraform bar:foo(a : T) where {T : Real}
            return a + self.index
        end
        terraform bar.sin(a : T, b ...) where {T : Real}
            return a + b._0 + 2
        end
        test bar.sin(2) == 4
        test mybar:foo(2) == 3
        test bar.sin(2,3) == 7
    end

    testset "concrete types" do
        terraform bar:foo(a : double, b : double)
            return a * b + self.index
        end
        terraform bar:foo(a : float, b : double)
            return a * b + self.index + 1
        end
        test mybar:foo(1., 2.) == 3.0
        test mybar:foo(float(1), 2.) == 4.0
    end

    testset "method dispatch parametric concrete types" do
        terraform bar:foo(a : T, b : T) where {T : Real}
            return a * b + self.index
        end
        terraform bar:foo(a : int, b : T) where {T : Real}
            return a * b + 1 + self.index
        end
        terraform bar:foo(a : T, b : int) where {T : Real}
            return a * b + 2 + self.index
        end
        terraform bar:foo(a : int, b : int)
            return a * b + 3 + self.index
        end
        --typical use-case in dispatch
        test mybar:foo(1., 2.) == 3.0
        test mybar:foo(1, 2.) == 4.0
        test mybar:foo(1., 2) == 5.0
        test mybar:foo(1, 2) == 6.0
    end

    testset "method dispatch parametric reference types" do
        terraform bar:foo(a : &T, b : &T) where {T : Real}
            return @a * @b + self.index
        end
        terraform bar:foo(a : &int, b : &T) where {T : Real}
            return @a * @b + 1 + self.index
        end
        terraform bar:foo(a : &T, b : &int) where {T : Real}
            return @a * @b + 2 + self.index
        end
        terraform bar:foo(a : &int, b : &int)
            return @a * @b + 3 + self.index
        end
        terracode
            var x = 1.
            var y = 2.
            var i = 1
            var j = 2
        end
        test mybar:foo(&x, &y) == 3.0
        test mybar:foo(&i, &y) == 4.0
        test mybar:foo(&x, &j) == 5.0
        test mybar:foo(&i, &j) == 6.0
    end

    testset "method dispatch mixed parametric concrete/reference types" do
        terraform bar:foo(a : T, b : &T) where {T : Real}
            return a + @b + self.index
        end
        terraform bar:foo(a : &T, b : T) where {T : Real}
            return @a + b + 1 + self.index
        end
        terracode
            var x = 1.
            var y = 2.
        end
        test mybar:foo(x, &y) == 4.0
        test mybar:foo(&x, y) == 5.0
    end

    testset "nested reference types" do
        terraform bar:foo(a : &&T, b : &&T) where {T : Real}
            return @(@a) * @(@b) + self.index
        end
        terraform bar:foo(a : &&int, b : &&T) where {T : Real}
            return @(@a) * @(@b) + 1 + self.index
        end
        terraform bar:foo(a : &&T, b : &&int) where {T : Real}
            return @(@a) * @(@b) + 2 + self.index
        end
        terraform bar:foo(a : &&int, b : &&int)
            return @(@a) * @(@b) + 3 + self.index
        end
        terracode
            var x, y, i, j = 1., 2., 1, 2
            var rx, ry, ri, rj = &x, &y, &i, &j
        end
        test mybar:foo(&rx, &ry) == 3.0
        test mybar:foo(&ri, &ry) == 4.0
        test mybar:foo(&rx, &rj) == 5.0
        test mybar:foo(&ri, &rj) == 6.0
    end

    testset "nearly ambiguous calls" do
        --foo<T>
        terraform bar:foo(a : T, b : T) where {T : Float}
            return a * b + self.index
        end
        --foo<T1,T2>
        terraform bar:foo(a : T1, b : T2) where {T1 : Float, T2 : Float}
            return a * b + 1 + self.index
        end
        --both of the following calls satisfy the concept dispatch 
        --checks on all arguments. 
        test mybar:foo(1., 2.) == 3.0 --calling foo<T> and foo<T1,T2> are both valid. However,
        --foo<T> is considered to be more specialized, so we pick this method.
        test mybar:foo(float(1), 2.) == 4.0 --calling foo<T> would lead to a cast, which is 
        --not allowed, so we call foo<T1,T2>
    end

    testset "nested namespaces" do
        local ns = {}
        ns.Float = Float
        ns.sin = {}
        ns.sin.Float = Float
        ns.sin.bar = bar
        --foo<T>
        terraform ns.sin.bar:foo(a : T, b : T) where {T : ns.sin.Float}
            return a * b + self.index
        end
        --foo<T1,T2>
        terraform ns.sin.bar:foo(a : T1, b : T2) where {T1 : ns.Float, T2 : ns.sin.Float}
            return a * b + 1  + self.index
        end
        terracode
            var mybar = ns.sin.bar{1}
        end
        test mybar:foo(1., 2.) == 3.0
        test mybar:foo(float(1), 2.) == 4.0 
    end

    testset "use in parametric class definition" do
        local Bar = function(T)
            local struct bar{
                x : T
            }
            terraform bar:eval(x : S) where {S}
                return x + self.x
            end
            return bar
        end
        local bbar = Bar(double)
        terracode
            var mybar = bbar{1}
        end
        test mybar:eval(1.) == 2
    end

end
