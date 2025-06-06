-- SPDX-FileCopyrightText: 2024 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2024 Torsten Keßler <t.kessler@posteo.de>
-- SPDX-FileCopyrightText: 2025 René Hiemstra <rrhiemstar@gmail.com>
-- SPDX-FileCopyrightText: 2025 Torsten Keßler <t.kessler@posteo.de>
--
-- SPDX-License-Identifier: MIT

require "terralibext"

local C = terralib.includecstring[[
    #include <stdio.h>
    #include <string.h>
]]

local base = require("base")
local interface = require("interface")
local range = require("range")
local err = require("assert")

local size_t = uint64
local u8 = uint8


local function Base(block, T, options)

    assert(
        options and type(options.copyable)=="boolean",
        "Invalid option. Please provide {copyable = true / false}"
    )
    --is the type copyable? Default to false.
    local copyable = options.copyable or false

    --type traits
    block.isblock = true
    block.type = block
    block.traits.eltype = T
    block.elsize = T==opaque and 1 or sizeof(T)

    block.methods.getdataptr = terra(self : &block)
        return self.ptr
    end

    --block is empty, no resource and no allocator
    block.methods.isempty = terra(self : &block)
        return self.ptr==nil
    end

    --resource is borrowed, there is no allocator
    --this represents a view of the data
    block.methods.borrows_resource = terra(self : &block)
        return self.ptr~=nil and self.alloc.data==nil
    end

    --resource is owned, there is an allocator
    block.methods.owns_resource = terra(self : &block)
        return self.ptr~=nil and self.alloc.data~=nil
    end

    block.methods.size_in_bytes = terra(self : &block) : size_t
        return self.nbytes
    end

    if T==opaque then
        block.methods.size = terra(self : &block) : size_t
            return self.nbytes
        end
    else
        block.methods.size = terra(self : &block) : size_t
            return self.nbytes / [block.elsize]
        end
    end

    --initialize to empty block
    block.methods.__init = terra(self : &block)
        self.ptr = nil
        self.nbytes = 0
        self.alloc.data = nil
        self.alloc.ftab = nil
    end

    --exact clone of the block
    block.methods.clone = terra(self : &block)
        var newblk : block --allocate memory for exact clone
        if not self:isempty() then
            self.alloc:__allocators_best_friend(&newblk, [ block.elsize ], self:size())
            if not newblk:isempty() then
                C.memcpy(newblk.ptr, self.ptr, self:size_in_bytes())
            end
        end
        return newblk
    end

end

--abstraction of a memory block without any type information.
local struct block

local __Allocator = interface.newinterface("__Allocator")
terra __Allocator:__allocators_best_friend(blk: &block, elsize: size_t, counter: size_t) end
__Allocator:complete()

struct block{
    ptr : &opaque
    nbytes : size_t
    alloc : __Allocator
}

function block.metamethods.__typename(self)
    return "block"
end

--add base functionality
base.AbstractBase(block)
Base(block, opaque, {copyable=false})

--__dtor for opaque memory block
terra block.methods.__dtor(self : &block)
    if self:borrows_resource() then
        self:__init()
    elseif self:owns_resource() then
        self.alloc:__allocators_best_friend(self, 0, 0)
    end
end

--add raii move method
terralib.ext.addmissing.__move(block)
block:complete()


--abstraction of a memory block with type information.
local SmartBlock
SmartBlock = terralib.memoize(function(T, options)

    --check optional input
    assert(options and type(options.copyable)=="boolean",
        "Invalid option. Please provide {copyable = true / false}"
    )

    --is the type copyable? Default to false.
    local copyable = options.copyable or false

    local struct block{
        ptr : &T
        nbytes : size_t
        alloc : __Allocator
    }

    function block.metamethods.__typename(self)
        return ("SmartBlock(%s)"):format(tostring(T))
    end

    base.AbstractBase(block)

    -- Cast block from one type to another
    function block.metamethods.__cast(from, to, exp)
        local function passbyvalue(to, from)
            if from:ispointertostruct() and to:ispointertostruct() then
                return false, to.type, from.type
            end
            return true, to, from
        end
        --process types
        local byvalue, to, from = passbyvalue(to, from)        
        --exit early if types do not match
        if not to.isblock or not from.isblock then
            error("Arguments to cast need to be of generic type SmartBlock.")
        end
        --perform cast
        if byvalue then
            --case when to.eltype is a managed type (implements a '__dtor')
            if terralib.ext.ismanaged(to.traits.eltype) then
                return quote
                    --we get a handle to the object, which means we get an lvalue that 
                    --does not own the resource, so it's '__dtor' will not be called
                    var tmp = __handle__(exp)
                    --debug check if sizes are compatible, that is, is the
                    --remainder zero after integer division
                    --err.assert(tmp:size_in_bytes() % [to.elsize]  == 0)
                    --loop over all elements of blk and initialize their entries 
                    var size = tmp:size_in_bytes() / [to.elsize]
                    var ptr = [&to.traits.eltype](tmp.ptr)
                    for i = 0, size do
                        ptr:__init()
                        ptr = ptr + 1
                    end
                in
                    [to.type]{[&to.traits.eltype](tmp.ptr), tmp.nbytes, tmp.alloc}
                end
            --simple case when to.eltype is not managed
            else
                return quote
                    var tmp = __handle__(exp)
                    --debug check if sizes are compatible, that is, is the
                    --remainder zero after integer division
                    --err.assert(tmp:size_in_bytes() % [to.elsize]  == 0)
                in
                    [to.type]{[&to.traits.eltype](tmp.ptr), tmp.nbytes, tmp.alloc}
                end
            end
        else
            --passing by reference
            return quote
                --store the reference such that we can access it.
                var blk = exp
                err.assert(blk:size_in_bytes() % [to.elsize]  == 0)
            in
                [&to.type](blk)
            end
        end
    end --__cast

    --declaring __dtor for use in implementation below
    terra block.methods.__dtor :: {&block} -> {}

    function block.metamethods.__staticinitialize(self)

        --add base functionality
        Base(block, T, options)

        --setters and getters
        block.methods.get = terra(self : &block, i : size_t)
            err.assert(i < self:size())
            return self.ptr[i]
        end

        block.methods.set = terra(self : &block, i : size_t, v : T)
            err.assert(i < self:size())
            self.ptr[i] = v
        end

        block.metamethods.__apply = macro(function(self, i)
            return quote
                err.assert(i < self:size())
            in
                self.ptr[i]
            end
        end)

        block.staticmethods.frombuffer = terra(size: size_t, ptr: &T)
            var nbytes = size * sizeof(T)
            var b: block
            b.ptr = ptr
            b.nbytes = nbytes
            b.alloc.data = nil
            b.alloc.ftab = nil
            return b
        end

        --iterator - behaves like a pointer and can be passed
        --around like a value, convenient for use in ranges.
        local struct iterator{
            parent : &block
            ptr : &T
        }

        terra block:getiterator()
            return iterator{self, self.ptr}
        end

        terra iterator:getvalue()
            return @self.ptr
        end

        terra iterator:next()
            self.ptr = self.ptr + 1
        end

        terra iterator:isvalid()
            return (self.ptr - self.parent.ptr) * [block.elsize] < self.parent.nbytes
        end
        
        block.iterator = iterator
        range.Base(block, iterator)

        terra block:reallocate(size : size_t)
            self.alloc:__allocators_best_friend(self, sizeof(T), size)
        end
        
        --implementation __dtor
        --ToDo: change recursion to a loop
        terra block.methods.__dtor(self : &block)
            --insert metamethods.__dtor if defined, which is used to introduce
            --side effects (e.g. counting number of calls for the purpose of testing)
            escape
                if block.metamethods and block.metamethods.__dtor then
                    emit quote
                        [block.metamethods.__dtor](self)
                    end
                end
            end
            --return if block is empty
            if self:isempty() then
                return
            end
            --initialize block and return if block borrows a resource
            if self:borrows_resource() then
                self:__init()
                return
            end
            --first destroy other memory block resources pointed to by self.ptr
            --ToDo: change recursion into a loop
            escape
                if terralib.ext.ismanaged(T) then
                    emit quote
                        var ptr = self.ptr
                        for i = 0, self:size() do
                            ptr:__dtor()
                            ptr = ptr + 1
                        end
                    end
                end
            end
            --free current memory block resources
            self.alloc:__allocators_best_friend(self, 0, 0)
        end

        --conditional compilation of a copy-method
        if copyable then
            block.methods.__copy = terra(from : &block, to : &block)
                --to:__dtor() is injected here by the compiler
                @to = from:clone()
            end
        end

        --add raii move method
        terralib.ext.addmissing.__move(block)

    end --__staticinitialize

	return block
end)


return {
    block = block,
    SmartBlock = SmartBlock
}
