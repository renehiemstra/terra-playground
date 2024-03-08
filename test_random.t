local random = require("random")
local C = terralib.includec("stdio.h")

terra main()
	escape
		for name, gen in pairs(random) do
			local rand = gen(double)
			emit quote
				var rng = [rand].new()
				var n: int64 = 2000001
				var mean: double = 0
				for i: int64 = 0, n do
					var u = rng:rand_uniform()
					mean = i * mean + u
					mean = mean / (i + 1)
				end
				C.printf("%s %u %g\n", name, n, mean)
			end
		end
	end
end

main()

