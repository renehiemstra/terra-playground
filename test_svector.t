import "terratest/terratest"

local tmath = require("mathfuns")
local svector = require("svector")


testenv "Static vector - construction" do

    for _,T in ipairs{int, int64, float, double} do
        for N=2,3 do

            local svec = svector.StaticVector(T,N)   
         
            testset(N,T) "new, size, get, set" do
                terracode
                    var v = svec.new()
                    for i = 0, N do              
                        v:set(i, i+1)
                    end                     
                end
                test v:size()==N
                for i=0,N-1 do              
                    test v:get(i) == T(i+1)
                end 
            end
          
            testset(N,T) "zeros" do                       
                terracode                                  
                    var v = svec.zeros()
                end
                test v:size()==N
                for i = 0, N-1 do              
                    test v:get(i) == 0
                end 
            end 
        
            testset(N,T) "ones" do                       
                terracode                                  
                    var v = svec.ones()
                end 
                test v:size()==N
                for i = 0, N-1 do              
                    test v:get(i) == T(1)
                end 
            end 

            testset(N,T) "fill" do                       
                terracode                                  
                    var v = svec.fill(T(3))
                end 
                test v:size()==N
                for i = 0, N-1 do              
                    test v:get(i) == T(3)
                end 
            end
            
        end --N

        testset "from (N=2)" do
            local svec = svector.StaticVector(T, 2)
            terracode
                var u = svec.from(1, 2)
                var v = svec.from(1, 2)
                var w = svec.from(0, 2)
            end
            test u:size() == 2
            test u:get(0) == 1
            test u:get(1) == 2
            test u == v
            test u ~= w
        end

        testset "from (N=3)" do
            local svec = svector.StaticVector(T, 3)
            terracode
                var u = svec.from(1, 2, 3)
                var v = svec.from(1, 2, 3)
                var w = svec.from(0, 2, 3)
            end
            test u:size() == 3
            test u:get(0) == 1
            test u:get(1) == 2
            test u:get(2) == 3
            test u == v
            test u ~= w
        end

    end --T

end --testenv


testenv "Static vector - arithmatic" do

    for _,T in ipairs{int, int64, float, double} do
        for N = 2,3 do
    
            local svec = svector.StaticVector(T,N)   

            terracode 
                var a = svec.fill(2)
                var b = svec.new()
                for i = 0, N do
                    b(i) = i+1
                end
            end

            testset(N,T) "__add" do
                terracode
                    var c = a + b
                    var d = b + 1
                    var e = 1 + b
                end
                for i=0, N-1 do          
                    test c(i) == 2 + i + 1
                    test d(i) == i + 2
                    test e(i) == i + 2
                end 
            end

            testset(N,T) "__sub" do
                terracode
                    var c = a - b  
                    var d = b - 1
                    var e = 1 - b
                end
                for i=0, N-1 do          
                    test c(i) == 2 - (i + 1)
                    test d(i) == i
                    test e(i) == -i
                end 
            end

            testset(N,T) "__mul" do
                terracode
                    var c = a * b
                    var d = b * 2
                    var e = 2 * b
                end
                for i=0, N-1 do          
                    test c(i) == 2*(i+1)
                    test d(i) == 2*(i+1)
                    test e(i) == 2*(i+1)
                end 
            end

            if not T:isintegral() then
                -- / doesnt work well on vector types that are not integral.
                testset(N,T) "__div" do
                    terracode
                        var c = a / b
                        var d = b / 2
                        var e = 2 / b
                    end
                    for i=0, N-1 do
                        test c(i) == T(2) / (T(i)+T(1))
                        test d(i) == (T(i)+T(1)) / T(2)
                        test e(i) == T(2) / (T(i)+T(1))
                    end 
                end
            end

            testset(N,T) "__sum" do
                for i=0, N-1 do          
                    test a:sum() == 2 * N
                end 
                if N % 2 == 0 then
                    for i=0, N-1 do          
                        test b:sum() == (N+1) * (N / 2)
                    end 
                else
                    for i=0, N-1 do          
                        test b:sum() == (N+1) * (N / 2) + (N+1)/2
                    end
                end
            end

        end --for N
    end --for _,T

end

testenv "Static vector - reductions" do

    local function calcnorm2(n)
        local s = 0
        for i = 1, n do
            s = s + i^2
        end
        return s
    end

    for _,T in ipairs{int, int64, float, double} do

        for N = 2,4 do
    
            local svec = svector.StaticVector(T,N)   

            terracode 
                var a = svec.fill(2)
                var b = svec.new()
                for i = 0, N do
                    b(i) = i+1
                end
            end

            testset(N,T) "sum" do
                test a:sum() == 2 * N
                if N % 2 == 0 then        
                    test b:sum() == (N+1) * (N / 2)
                else       
                    test b:sum() == (N+1) * (N / 2) + (N+1)/2
                end
            end

            testset(N,T) "norm2" do
                test a:norm2() == 4 * N
                test b:norm2() == [ calcnorm2(N) ]
            end

            testset(N,T) "norm" do
                terracode
                    var a_norm = a:norm()
                end
                test tmath.isapprox(a:norm() * a:norm(), a:norm2(), 1e-14)
                test tmath.isapprox(b:norm() * b:norm(), b:norm2(), 1e-14)
            end

        end --for N
    end --for _,T

end
