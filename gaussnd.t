local io = terralib.includec('stdio.h')
local alloc = require('alloc')
local gauss = require("gauss")
local rn = require("range")



local Allocator = alloc.Allocator
local DefaultAllocator =  alloc.DefaultAllocator()

local T = double
local N = 10

local terra main()
    var alloc : DefaultAllocator
    var x1, w1 = gauss.legendre(&alloc, N)
    var x2, w2 = gauss.legendre(&alloc, N)
    for u in rn.product(x1, x2) do
        io.printf("x = (%0.2f, %0.2f)\n", u._0, u._1)
    end
end
main()