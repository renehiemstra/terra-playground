-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
-- SPDX-FileCopyrightText: 2025 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2025 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

import "terraform"

local base = require("base")
local concepts = require("concepts")
local template = require("template")
local lambda = require("lambda")
local tmath = require("tmath")
local nfloat = require("nfloat")
local err = require("assert")
local parametrized = require("parametrized")

local size_t = uint64

--get the terra-type of a pointer or type
local gettype = function(t)
    assert(terralib.types.istype(t) and "Not a terra type")
    if t:ispointer() then
        return t.type
    else
        return t
    end
end

--given a terra value or reference to a value, get its value
local byvalue = function(t)
    local typ = t.type or t.tree.type or error("Not a terra type.")
    if typ:ispointer() then
        return `@t
    else
        return `t
    end
end

--return the value-type of an iterator type
local function getvalue_t(iterator_t)
    return iterator_t.methods.getvalue.type.returntype
end

local IteratorBase = function(Iterator)
    --type trait value type
    Iterator.value_t = getvalue_t(Iterator)
end

--an iterator implements the following macros:
--  methods.getfirst :: (self) -> (state, value)
--  methods.getnext :: (self, state) -> (value)
--  methods.isvalid :: (self, state, value) -> (true/false)
--the following base class then overloads the '>>' operator
--and adds the '__for' metamethod, and adds a 'collect' 
--method that collects all elements in the range in a container
--that satsifies the 'Stacker(T)' interface
local RangeBase = function(Range, iterator_t)

    --set base functionality for iterators
    IteratorBase(iterator_t)

    --set the value type and iterator type of the range
    Range.isrange = true
    Range.iterator_t = iterator_t
    Range.value_t = iterator_t.value_t

    --overloading '>>' operator
    Range.metamethods.__rshift = macro(function(self, adapter)
        local self_type = self.tree.type
        local adapter_type = adapter.tree.type
        if self_type:isstruct() and self_type.metamethods.__for 
            and adapter_type:isstruct() 
        then
            local Adapter = adapter_type.generator
            local A = Adapter(self_type, adapter_type)
            -- HACK Trigger init() manually
            -- Normally, we would do return `A {self, adapter}
            -- However, the current RAII implementation does not cover all
            -- possible cases for an init() call yet.
            -- Below code first declares an element of type A and then
            -- fills the entries of A by the corresponding value. The situation
            -- is complicated by the fact that each adaptor has different
            -- entries names. Hence, we need to do a little bit of meta
            -- programming to extract the name of the struct field.
            assert(#A.entries == 2)
            local newrange = symbol(A)
            return (
                quote 
                    escape
                        emit quote var [newrange] end
                        for i, v in ipairs{self, adapter} do
                            local name = A.entries[i].field
                            emit quote [newrange].[name] = v end
                        end
                    end
                in
                    [newrange]
                end
            )
        end
    end)

    --__for is generated for iterators
    Range.metamethods.__for = function(self, body)
        return quote
            var iter = self:getiterator()
            while iter:isvalid() do             --while not at the end
                var value = iter:getvalue()     --get value
                [body(byvalue(value))]          --run body of loop
                iter:next()                     --increment state
            end
        end
    end

    local T = Range.value_t
    local Stack = concepts.Stack(T)
    local DStack = concepts.DStack(T)
    
    --containers that only implement the setter interface are using 'set'. Sufficient
    --space needs to be allocated before
    terraform Range:collect(container : &S) where {S : Stack}
        var i = 0
        for v in self do
            container:set(i, v)
            i = i + 1
        end
    end

    --dynamic stacks have a push
    terraform Range:pushall(container : &S) where {S : DStack}
        for v in self do
            container:push(v)
        end
    end

end

local floor
terraform floor(v : T) where {T : concepts.Integer}
    return size_t(v)
end

terraform floor(v : T) where {T : concepts.Float}
    return [size_t](tmath.floor(v))
end

terraform floor(v : T) where {T : concepts.NFloat}
    return [size_t](v:truncatetodouble())
end

local truncate
terraform truncate(v : T) where {T}
    return [size_t](v)
end

terraform truncate(v : T) where {T : concepts.NFloat}
    return [size_t](v:truncatetodouble())
end


local unitrange = parametrized.type(function(T)

    local struct range{
        a : T
        b : T
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(range)

    local new = terra(a : T, b : T, include_last : bool)
        err.assert((b-a) > 0)
        var size = floor(b-a) + [int](include_last)
        return range{a, a + size}
    end

    range.staticmethods.new = terralib.overloadedfunction("new",{
        new,
        terra(a : T, b : T) return new(a, b, false) end
    })

    terra range:length()
        return self.b - self.a
    end

    range.metamethods.__apply = terra(self : &range, i : size_t)
        err.assert(i < self:length())
        return self.a + i
    end

    local struct iterator{
        parent : &range
        state : T
    }

    terra iterator:next()
        self.state = self.state + 1
    end

    terra iterator:getvalue()
        return self.state
    end

    terra iterator:isvalid()
        return self.state < self.parent.b
    end

    terra range:getiterator()
        return iterator{self, self.a}
    end

    --add metamethods
    RangeBase(range, iterator)

    return range
end)

local steprange = parametrized.type(function(T)

    local struct range{
        a : T
        b : T
        step : T
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(range)

    local new = terra(a : T, b : T, step : T, include_last : bool)
        err.assert(((b-a) >= 0 and step > 0) or ((b-a) <= 0 and step < 0))
        b = terralib.select(b > a, b + [int](include_last), b - [int](include_last))
        b = b + (b - a) % step
        return range{a, b, step}
    end

    range.staticmethods.new = terralib.overloadedfunction("new",{
        new,
        terra(a : T, b : T, step : T) return new(a, b, step, false) end
    })
    
    terra range:length() : size_t
        return truncate((self.b-self.a) / self.step)
    end

    range.metamethods.__apply = terra(self : &range, i : size_t)
        err.assert(i < self:length())
        return self.a + i * self.step
    end

    local struct iterator{
        parent : &range
        state : T
    }

    terra iterator:next()
        self.state = self.state + self.parent.step
    end

    terra iterator:getvalue()
        return self.state
    end

    terra iterator:isvalid()
        return terralib.select(self.parent.step>0, self.state < self.parent.b, self.state > self.parent.b)
    end

    terra range:getiterator()
        return iterator{self, self.a}
    end

    --add metamethods
    RangeBase(range, iterator)

    return range
end)

local infunitrange = parametrized.type(function(T)

    local struct range{
        a : T
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(range)

    range.staticmethods.new = terra(a : T)
        return range{a}
    end

    range.metamethods.__apply = terra(self : &range, i : size_t)
        return self.a + i
    end

    local struct iterator{
        parent : &range
        state : T
    }

    terra iterator:next()
        self.state = self.state + 1
    end

    terra iterator:getvalue()
        return self.state
    end

    terra iterator:isvalid()
        return true
    end

    terra range:getiterator()
        return iterator{self, self.a}
    end

    --add metamethods
    RangeBase(range, iterator)

    return range
end)

local infsteprange = parametrized.type(function(T)

    local struct range{
        a : T
        step : T
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(range)

    range.staticmethods.new = terra(a : T, step : T)
        return range{a, step}
    end
    
    range.metamethods.__apply = terra(self : &range, i : size_t)
        return self.a + i * self.step
    end

    local struct iterator{
        parent : &range
        state : T
    }

    terra iterator:next()
        self.state = self.state + self.parent.step
    end

    terra iterator:getvalue()
        return self.state
    end

    terra iterator:isvalid()
        return true
    end

    terra range:getiterator()
        return iterator{self, self.a}
    end

    --add metamethods
    RangeBase(range, iterator)

    return range
end)

local Unitrange = parametrized.type(function(T, sentinal)
    local sentinal = sentinal or "bounded"
    if sentinal == "bounded" then
        return unitrange(T)
    elseif sentinal == "infinite" then
        return infunitrange(T)
    else
        error("ArgumentError: second (optional) argument should be 'bounded' or 'infinite'.")
    end
end)

local Steprange = parametrized.type(function(T, sentinal)
    local sentinal = sentinal or "bounded"
    if sentinal == "bounded" then
        return steprange(T)
    elseif sentinal == "infinite" then
        return infsteprange(T)
    else
        error("ArgumentError: second (optional) argument should be 'bounded' or 'infinite'.")
    end
end)

local TransformedRange = function(Range, Function)

    local struct transform{
        range : Range
        f : Function
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(transform)

    local iterator_t = Range.iterator_t
    local T = Range.value_t

    local eval = macro(function(self, value)
        if T.convertible=="tuple" then --we always unpack tuples
            return quote 
                var v = value
            in
                self.f(unpacktuple(v))
            end
        else
            return `self.f(value)
        end
    end)
    
    local struct iterator{
        adapter : &transform
        state : iterator_t
    }

    terra iterator:next()
        self.state:next()
    end

    terra iterator:getvalue()
        return eval(self.adapter, self.state:getvalue())
    end

    terra iterator:isvalid()
        return self.state:isvalid()
    end

    terra transform:getiterator()
        return iterator{self, self.range:getiterator()}
    end

    --add metamethods
    RangeBase(transform, iterator)

    return transform
end

local FilteredRange = function(Range, Function)

    local struct filter{
        range : Range
        predicate : Function
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(filter)

    local iterator_t = Range.iterator_t
    local T = Range.value_t

    --evaluate predicate
    local pred = macro(function(self, value)
        if T.convertible=="tuple" then --we always unpack tuples
            return quote 
                var v = value
            in
                self.predicate(unpacktuple(v))
            end
        else
            return `self.predicate(value)
        end
    end)

    local struct iterator{
        adapter : &filter
        state : iterator_t
    }

    terra filter:getiterator()
        var state = self.range:getiterator()
        while (state:isvalid() and pred(self, state:getvalue())==false) do
            state:next()
        end
        return iterator{self, state}
    end

    terra iterator:next()
        repeat
            self.state:next()
        until pred(self.adapter, self.state:getvalue()) or not self.state:isvalid()
    end

    terra iterator:getvalue()
        return self.state:getvalue()
    end

    terra iterator:isvalid()
        return self.state:isvalid()
    end

    --add metamethods
    RangeBase(filter, iterator)

    return filter
end

local TakeRange = function(Range)

    local struct take{
        range : Range
        count : int64
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(take)

    local iterator_t = Range.iterator_t

    local struct iterator{
        adapter : &take
        state : iterator_t
        count : int64
    }

    terra take:getiterator()
        return iterator{self, self.range:getiterator(), self.count}
    end

    terra iterator:next()
        self.state:next()
        self.count = self.count - 1
    end

    terra iterator:getvalue()
        return self.state:getvalue()
    end

    terra iterator:isvalid()
        return self.state:isvalid() and self.count > 0
    end

    --add metamethods
    RangeBase(take, iterator)

    return take
end

local DropRange = function(Range)

    local struct drop{
        range : Range
        count : int64
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(drop)

    local iterator_t = Range.iterator_t

    local struct iterator{
        adapter : &drop
        state : iterator_t
    }

    terra drop:getiterator()
        var state = self.range:getiterator()
        var count = 0
        while state:isvalid() and count < self.count do
            state:next()
            count = count + 1
        end
        return iterator{self, state}
    end

    terra iterator:next()
        self.state:next()
    end

    terra iterator:getvalue()
        return self.state:getvalue()
    end

    terra iterator:isvalid()
        return self.state:isvalid()
    end

    --add metamethods
    RangeBase(drop, iterator)

    return drop
end

local TakeWhileRange = function(Range, Function)

    local struct takewhile{
        range : Range
        predicate : Function
    }
    --add methods, staticmethods and templates and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(takewhile)

    local iterator_t = Range.iterator_t
    local T = Range.value_t

    --evaluate predicate
    local pred = macro(function(self, value)
        if T.convertible=="tuple" then --we always unpack tuples
            return quote 
                var v = value
            in
                self.predicate(unpacktuple(v))
            end
        else
            return `self.predicate(value)
        end
    end)

    local struct iterator{
        adapter : &takewhile
        state : iterator_t
    }

    terra takewhile:getiterator()
        return iterator{self, self.range:getiterator()}
    end

    terra iterator:next()
        self.state:next()
    end

    terra iterator:getvalue()
        return self.state:getvalue()
    end

    terra iterator:isvalid()
        return self.state:isvalid() and pred(self.adapter,self.state:getvalue())==true
    end

    --add metamethods
    RangeBase(takewhile, iterator)

    return takewhile
end

local DropWhileRange = function(Range, Function)

    local struct dropwhile{
        range : Range
        predicate : Function
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(dropwhile)

    local iterator_t = Range.iterator_t
    local T = Range.value_t

    --evaluate predicate
    local pred = macro(function(self, value)
        if T.convertible=="tuple" then --we always unpack tuples
            return quote 
                var v = value
            in
                self.predicate(unpacktuple(v))
            end
        else
            return `self.predicate(value)
        end
    end)

    local struct iterator{
        adapter : &dropwhile
        state : iterator_t
    }

    terra dropwhile:getiterator()
        var state = self.range:getiterator()
        while state:isvalid() and pred(self,state:getvalue()) do
            state:next()
        end
        return iterator{self, state}
    end

    terra iterator:next()
        self.state:next()
    end

    terra iterator:getvalue()
        return self.state:getvalue()
    end

    terra iterator:isvalid()
        return self.state:isvalid()
    end

    --add metamethods
    RangeBase(dropwhile, iterator)

    return dropwhile
end

--factory function for range adapters that carry a lambda
local adapter_lambda_factory = function(Adapter)
    local factory = macro(
        function(fun, capture)
            --get the captured variables
            local p = lambda.makelambda(fun, capture or `{})
            --set the generator (FilteredRange or TransformedRange, etc)
            local lambda_t = p:gettype()
            lambda_t.generator = Adapter
            --create and return lambda object by value
            return `p
        end)
    return factory
end

--factory function for range adapters that carry a view
local adapter_view_factory = function(Adapter)
    local factory = macro(
        function(n)
            --wrapper struct
            local struct view{
                size : int64
            }
            view.generator = Adapter
            --enable casting to an integer from view
            view.metamethods.__cast = function(from, to, exp)
                if from:isstruct() and to:isintegral() then
                    return quote
                        var x = exp
                    in
                        [int64](x.size)
                    end
                end
            end
            --create and return wrapper object by value
            return `view {n}
        end)
    return factory
end

--factory function for range adapters that don't cary state
local adapter_simple_factory = function(Adapter)
    local factory = macro(function()
        --wrapper struct
        local struct simple{
        }
        simple.generator = Adapter
        --create and return simple object by value
        return `simple {}
    end)
    return factory
end

local function getunderlyingtype(t)
    if not terralib.types.istype(t) then
        t = t:gettype()
    end
    if t:ispointer() then
        return t.type
    else
        return t
    end
end

--factory function for range combiners
local combiner_factory = function(Combiner)
    local combiner = macro(function(...)
        --take all arguments
        local args = terralib.newlist{...}
        local N = #args
        --filter between ranges and options
        local ranges = args:filter(function(v) return getunderlyingtype(v).isrange end)
        local options
        if not args[N]:gettype().isrange then
            options = args[N]:asvalue()
        end
        --get range types
        local range_types = ranges:map(function(v) return v:gettype() end)
        --construct the combirange type and instantiate 
        --terra obj
        local combirange = Combiner(range_types, options)
        return quote
            --var range = combirange{[ranges]}
            var range : combirange
            --HACK: needed for now because initializers don't work correctly
            --with raii
            escape
                for i,e in ipairs(combirange:getentries()) do
                    emit quote
                        range.[e.field] = [ ranges[i] ]
                    end
                end
            end
        in
            range
        end
    end)
    return combiner
end

local newcombiner = function(Ranges, name)
    --create struct
    local combiner = terralib.types.newstruct(name)
    --add entries
    for i,Range in ipairs(Ranges) do
		combiner.entries:insert({field = "_"..tostring(i-1), type = Range})
	end
    combiner:setconvertible("tuple")
    --complete struct type
	combiner:complete()

    return combiner
end

local Enumerator = function(Ranges)

    --check that a range-for is implemented
    assert(#Ranges==1)
    local Range = Ranges[1]
    local byreference = Range:ispointer()
    if not (byreference and Range.type.metamethods.__for or Range.metamethods.__for) then
        error("Terra type does not implement the range interface.")
    end

    local struct enumerator{
        range : Range
    }

    enumerator.metamethods.__for = function(self,body)
        return quote
            var iter = self
            var i = 0
            for v in [byvalue(`iter.range)] do
                [body(i,v)]
                i = i + 1
            end
        end
    end

    return enumerator
end

local JoinRange = function(Ranges)

    local joiner = newcombiner(Ranges, "joiner")
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(joiner)
    local D = #Ranges

    --get value, range and iterator types
    local T = gettype(Ranges[1]).value_t
    local state_t = terralib.newlist{}
    for i,rn in ipairs(Ranges) do
        assert(gettype(rn).value_t == T, "ArgumentError: the value type should be uniform.") --make sure the value type is uniform
        state_t:insert(gettype(rn).iterator_t)
    end
    local iterator_t = tuple(unpack(state_t))

    local struct iterator{
        range : &joiner
        state : iterator_t
        index : uint8
    }
    
    terra joiner:getiterator()
        var iter : iterator
        iter.range = self
        escape
            for k=0,D-1 do
                local s = "_"..tostring(k)
                emit quote
                    iter.state.[s] = self.[s]:getiterator()
                end
            end
        end
        iter.index = 0
        return iter
    end

    terra iterator:getvalue()
        escape
            for k = 0, D-1 do
                local s = "_" .. tostring(k)
                emit quote
                    if self.index == [k] then
                        return self.state.[s]:getvalue()
                    end
                end
            end
        end
    end

    terra iterator:next()
        escape
            for k = 0, D-1 do
                local s = "_" .. tostring(k)
                emit quote
                    if self.index == [k] then
                        --advance iterator k
                        self.state.[s]:next()
                    end
                end
            end
        end
    end

    terra iterator:isvalid()
        escape
            for k = 0, D-1 do
                local s = "_" .. tostring(k)
                emit quote
                    if self.index == [k] then
                        if not self.state.[s]:isvalid() then
                            --jump to next iterator
                            self.index = self.index + 1
                        end
                    end
                end
            end
        end
        return self.index ~= D
    end

    --add metamethods
    RangeBase(joiner, iterator, T)
    
    return joiner
end

local ZipRange = function(Ranges)
  
    local zipper = newcombiner(Ranges, "zip")
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(zipper)
    local D = #Ranges

    --get range types
    local state_t = terralib.newlist{}
    local value_t = terralib.newlist{}
    for i,rn in ipairs(Ranges) do
        state_t:insert(gettype(rn).iterator_t)
        value_t:insert(gettype(rn).value_t)
    end
    local iterator_t = tuple(unpack(state_t))
    local T = tuple(unpack(value_t))

    local struct iterator{
        range : &zipper
        state : iterator_t
    }
    
    terra zipper:getiterator()
        var iter : iterator
        iter.range = self
        escape
            for k=0,D-1 do
                local s = "_"..tostring(k)
                emit quote
                    iter.state.[s] = self.[s]:getiterator()
                end
            end
        end
        return iter
    end

    terra iterator:getvalue()
        var value : T
        escape
            for k=0, D-1 do
                local s = "_"..tostring(k)
                emit quote value.[s] = self.state.[s]:getvalue() end
            end
        end
        return value
    end

    terra iterator:next()
        escape
            for k=0, D-1 do
                local s = "_"..tostring(k)
                emit quote self.state.[s]:next() end
            end
        end
    end

    terra iterator:isvalid()
        escape
            --loop over each of the D ranges
            for k=0, D-1 do
                local s = "_"..tostring(k)
                emit quote
                    if not self.state.[s]:isvalid() then
                        return false
                    end
                end
            end
        end
        return true
    end
    
    --add metamethods
    RangeBase(zipper, iterator)

    return zipper
end

local ProductRange = function(Ranges, options)

    --perm is a sequence of numbers denoting the perm in which the
    --product iterator iterates.

    local product = newcombiner(Ranges, "product")
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(product)
    local D = #Ranges
    --default ordering is {D, D-1, ... , 1}
    --this default is chosen to support array indexing in row-major format.
    local perm = options and options.perm or Ranges:mapi(function(i,v) return D+1-i end)
    assert(type(perm) == "table" and #perm == D)

    --local type traits
    local state_t = terralib.newlist{}
    local value_t = terralib.newlist{}
    for i,rn in ipairs(Ranges) do
        state_t:insert(gettype(rn).iterator_t)
        value_t:insert(gettype(rn).value_t)
    end
    local iterator_t = tuple(unpack(state_t))
    local T = tuple(unpack(value_t))

    --iterator definition
    local struct iterator{
        range : &product 
        state : iterator_t
        value : T
    }

    terra product:getiterator()
        var iter : iterator
        iter.range = self
        escape
            for k=0, D-1 do
                local s = "_"..tostring(k)
                emit quote 
                    iter.state.[s] = self.[s]:getiterator() 
                    iter.value.[s] = iter.state.[s]:getvalue()
                end
            end
        end
        return iter
    end

    terra iterator:getvalue()
        return self.value
    end

    terra iterator:next()
        escape
            for k=1, D-1 do
                local s = "_"..tostring(perm[k]-1)
                emit quote
                    --increase k
                    self.state.[s]:next()
                    if self.state.[s]:isvalid() then
                        self.value.[s] = self.state.[s]:getvalue()
                        return
                    end
                    --reset k
                    self.state.[s] = self.range.[s]:getiterator()
                    self.value.[s] = self.state.[s]:getvalue()
                end
            end
            --increase D-1
            local s = "_"..tostring(perm[D]-1)
            emit quote
                self.state.[s]:next()
                if self.state.[s]:isvalid() then
                    self.value.[s] = self.state.[s]:getvalue()
                end
            end
        end
    end

    terra iterator:isvalid()
        escape
            local s = "_"..tostring(perm[D]-1)
            emit quote
                return self.state.[s]:isvalid()
            end
        end
    end

    --add metamethods
    RangeBase(product, iterator)

    return product
end


local FoldLeft = function(Range, Function)
    
    --checking function input and outpiut arguments
    assert(Function.returntype, "ArgumentError: return type information is not available")
    assert(#Function.parameters == 2, "ArgumentError: not a binary function.")

    --pass by reference
    local passbyvalue = true
    if Function.returntype.convertible == "tuple" then
        assert(#Function.returntype.entries == 0)
        passbyvalue = false
    end

    --scalar datatype
    local T1 = Function.parameters[1]
    local T2 = Function.parameters[2]

    local struct foldl{
        range : Range
        f : Function
    }
    --add methods, staticmethods and templates tablet and template fallback mechanism 
    --allowing concepts-based function overloading at compile-time
    base.AbstractBase(foldl)

    if passbyvalue then -- pass-by-value
        foldl.methods.accumulatefrom = terra(self : &foldl, save : T1)
            for v in self.range do
                save = self.f(save, v)
            end
            return save
        end
    else -- pass-by-reference - T isa pointer
        foldl.methods.accumulatefrom = terra(self : &foldl, save : T1)
            for v in self.range do
                self.f(save, v)
            end
        end
    end

    return foldl
end

--generate user api macro's for adapters
local transform = adapter_lambda_factory(TransformedRange)
local filter = adapter_lambda_factory(FilteredRange)
local take = adapter_view_factory(TakeRange)
local drop = adapter_view_factory(DropRange)
local take_while = adapter_lambda_factory(TakeWhileRange)
local drop_while = adapter_lambda_factory(DropWhileRange)
--generate user api macro's for combi-ranges
local enumerate = combiner_factory(Enumerator)
local join = combiner_factory(JoinRange)
local product = combiner_factory(ProductRange)
local zip = combiner_factory(ZipRange)
--accumulators
local foldl = adapter_lambda_factory(FoldLeft)

--define reduction as a transform
local binaryoperation = {
    add = macro(function(x,y) return `x + y end),
    mul = macro(function(x,y) return `x * y end),
    div = macro(function(x,y) return `x / y end)
}

local reduce = macro(function(binaryop) 
    --reduction vararg template function
    local terraform tuplereduce(args ...)
        var res = args._0
        escape
            local n = #args.type.entries
            for i = 2, n do
                local s = "_" .. tostring(i-1)
                emit quote
                    res = binaryop(res, args.[s])
                end
            end
        end
        return res
    end
    return `transform(tuplereduce)
end)

local reverse = macro(function() 
    --reduction vararg template function
    local rev = macro(function(...)
        local args = terralib.newlist{...}
        local n = #args
        if n==1 and args[1]:gettype().convertible=="tuple" then
            --case of a tuple
            local reversedargs = terralib.newlist()
            for i = 1, n do
                local s = "_" .. tostring(n-i)
                reversedargs[i] = quote args.[s] end
            end
            return `{[reversedargs]}
        else
            return `{[args:rev()]}
        end
    end)
    --call transform to apply the above macro
    return `transform(rev)
end)

--math functions

terraform tmath.isapprox(A : &V, v : T, atol : S) where {V : concepts.Range, T : concepts.Number, S : concepts.Real}
    for a in A do
        if not tmath.isapprox(a, v, atol) then
            return false
        end
    end
    return true
end

terraform tmath.isapprox(A : &V, rn : &R, atol : S) where {V : concepts.Range, R : concepts.Range, S : concepts.Real}
    for t in zip(A, rn) do
        var a, v = t
        if not tmath.isapprox(a, v, atol) then
            return false
        end
    end
    return true
end

terraform tmath.isapprox(A : &V, rn : R, atol : S) where {V : concepts.Range, R : concepts.Range, S : concepts.Real}
    return tmath.isapprox(A, &rn, atol)
end


--export functionality for developing new ranges
local develop = {
    RangeBase = RangeBase,
    factory = { 
        combiner = combiner_factory,
        simple_adapter = adapter_simple_factory,
        view_adapter = adapter_view_factory,
        lambda_adapter = adapter_view_factory
    },
    newcombinerstruct = newcombiner,
}

--return module
return {
    include_last = true,
    Base = RangeBase,
    Unitrange = Unitrange,
    Steprange = Steprange,
    reduce = reduce,
    reverse = reverse,
    op = binaryoperation,
    transform = transform,    
    filter = filter,
    take = take,
    drop = drop, 
    take_while = take_while,
    drop_while = drop_while,
    enumerate = enumerate,
    join = join,
    product = product,
    zip = zip,
    foldl = foldl,
    develop = develop
}
