local io = terralib.includec("stdio.h")

local tmath = require("mathfuns")
local svector = require("svector")
local squad = require("singularquad")
local err = require("assert")

import "terratest/terratest"


testenv "interval" do

    local interval = squad.Interval(double)
    
    testset "new, vol, apply" do
        terracode
            var I = interval.new(1, 3)
        end
        test I.center==2
        test I.reach==1
        test I.origin == 1
        test I:vol()==2
    end

    testset "evaluation" do
        terracode
            var I = interval.new(1, 3)
            I.origin = 1.2
            var J = interval.new(1, 3)
            J.origin = 2.8
        end
        --test I(-0.1)==1 and I(0.9)==3
        test tmath.isapprox(I(-0.1), 1.0, 1e-15) 
        test tmath.isapprox(I(0.9), 3.0, 1e-15)
        test tmath.isapprox(J(-0.1), 3.0, 1e-15) 
        test tmath.isapprox(J(0.9), 1.0, 1e-15)
    end

    testset "barycentric coordinate - matching orientation" do
        terracode
            var I = interval.new(1, 3)
            I:setorigin(1)
        end
        test I:barycentriccoord(1)==0
        test I:barycentriccoord(2)==0.5
        test I:barycentriccoord(3)==1
    end

    testset "barycentric coordinate - reverse orientation" do
        terracode
            var I = interval.new(1, 3)
            I:setorigin(3)
        end
        test I:barycentriccoord(3)==0
        test I:barycentriccoord(2)==0.5
        test I:barycentriccoord(1)==1
    end

    testset "empty intersection" do
        terracode
            var I = interval.new(0, 1)
            var J = interval.new(2, 3)
            var Z = interval.intersect(I, J)
        end
        test Z.reach == -1
        test Z:isempty()
    end

    testset "nonempty intersection - interval" do
        terracode
            var I = interval.new(0, 2)
            var J = interval.new(1, 3)
            var Z = interval.intersect(I, J)
        end
        test Z==interval.new(1, 2)
        test Z.center==1.5
        test Z.reach==0.5
    end

    testset "nonempty intersection - point" do
        terracode
            var I = interval.new(0, 1)
            var J = interval.new(1, 2)
            var Z = interval.intersect(I, J)
        end
        test Z==interval.new(1, 1)
        test Z.center==1
        test Z.reach==0
    end

    testset "comparison" do
        terracode
            var I = interval.new(0, 2)
            var J = interval.new(0, 2)
            var Z = interval.new(1, 2)
        end
        test I == J
        test I ~= Z
    end
end

testenv "hypercube - N = 2" do

    local point = squad.Point(double, 2)
    local interval = squad.Interval(double)
    local hypercube = squad.Hypercube(double, 2)
    
    testset "new, dim, vol" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(0, 0))
            var C = hypercube.new(interval.new(0, 0), interval.new(0, 0))
        end
        test A:isempty()==false and B:isempty()==false and C:isempty()==false
        test A:dim()==2 and A:issingulardir(0)==false and A:issingulardir(1)==false
        test B:dim()==1 and B:issingulardir(0)==false and B:issingulardir(1)==true
        test C:dim()==0 and C:issingulardir(0)==true  and C:issingulardir(1)==true
        test A:vol()==4 and B:vol()==2 and C:vol()==1
    end

    testset "permutation" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5))
            var B = hypercube.new(interval.new(0, 2), interval.new(1, 1))
            var C = hypercube.new(interval.new(0, 0), interval.new(1, 2))
        end
        test A.perm[0] == 0 and A.perm[1] == 1
        test B.perm[0] == 0 and B.perm[1] == 1
        test C.perm[0] == 1 and C.perm[1] == 0
    end

    testset "square evaluation - origin==center" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5))
        end
        test A(0,0) == point.from(0,1)
        test A(0.5,0.5) == point.from(1,3)
        test A(1,1) == point.from(2,5)
    end

    testset "square evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5))
            A:setorigin(point.from(2,5))
        end
        test A(0,0) == point.from(2,5)
        test A(1,1) == point.from(0,1)
    end

    testset "segment evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(1, 3), interval.new(1, 1))
            A:setorigin(point.from(3,1))
        end
        test A(0) == point.from(3,1)
        test A(1) == point.from(1,1)
    end

    testset "square - barycentric coordinate" do
        terracode
            var A = hypercube.new(interval.new(1, 3), interval.new(2, 4))
            A:setorigin(point.from(1,2))
            var x = A:barycentriccoord(point.from(1,2))
        end
        test x._0==0 and x._1==0
    end

    testset "line - barycentric coordinate" do
        terracode
            var A = hypercube.new(interval.new(1, 3), interval.new(0, 0))
            A:setorigin(point.from(1,0))
            var x = A:barycentriccoord(point.from(2,0))
            var y = A:barycentriccoord(point.from(2,1))
        end
        test x._0==0.5 and x._1==0
        test y._0==0.5 and y._1==0 --projection
    end

    testset "empty intersection" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(3, 4), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()
    end

    testset "nonempty intersection - surface" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(1, 3), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(1, 2), interval.new(0, 2))
    end

    testset "nonempty intersection - line" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(2, 3), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(2, 2), interval.new(0, 2))
    end

    testset "nonempty intersection - point" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(2, 3), interval.new(2, 3))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(2, 2), interval.new(2, 2))
    end

    testset "product" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 0))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 2))
            var C = A * B
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(0, 2), interval.new(0, 2))
    end

    testset "division" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 0))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 2))
            var C = A * B
        end
        test (C / A) * A == C
        test (C / B) * B == C
    end
end

testenv "hypercube - N = 3" do

    local point = squad.Point(double, 3)
    local interval = squad.Interval(double)
    local hypercube = squad.Hypercube(double, 3)
    
    testset "new, dim, vol" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 0))
            var C = hypercube.new(interval.new(0, 2), interval.new(0, 0), interval.new(0, 0))
            var D = hypercube.new(interval.new(0, 0), interval.new(0, 0), interval.new(0, 0))
        end
        test A:isempty()==false and B:isempty()==false and C:isempty()==false and D:isempty()==false
        test A:dim()==3 and A:issingulardir(0)==false and A:issingulardir(1)==false and A:issingulardir(2)==false
        test B:dim()==2 and B:issingulardir(0)==false and B:issingulardir(1)==false and B:issingulardir(2)==true
        test C:dim()==1 and C:issingulardir(0)==false and C:issingulardir(1)==true  and C:issingulardir(2)==true
        test D:dim()==0 and D:issingulardir(0)==true  and D:issingulardir(1)==true  and D:issingulardir(2)==true
        test A:vol()==8 and B:vol()==4 and C:vol()==2 and D:vol()==1
    end

    testset "permutation" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 2), interval.new(0, 2))
            var C = hypercube.new(interval.new(0, 0), interval.new(0, 1), interval.new(0, 0))
            var D = hypercube.new(interval.new(0, 0), interval.new(0, 0), interval.new(0, 0))
        end
        test A.perm[0] == 0 and A.perm[1] == 1 and A.perm[2] == 2
        test B.perm[0] == 1 and B.perm[1] == 2 and B.perm[2] == 0
        test C.perm[0] == 1 and C.perm[1] == 2 and C.perm[2] == 0
        test D.perm[0] == 2 and D.perm[1] == 1 and D.perm[2] == 0
    end

    testset "barycentric coordinate - cube" do
        terracode
            var A = hypercube.new(interval.new(1, 3), interval.new(2, 4), interval.new(0, 1))
            A:setorigin(point.from(3,4,1))
            var x = A:barycentriccoord(point.from(3,4,1))
            var y = A:barycentriccoord(point.from(1,2,0))
        end
        test x._0==0 and x._1==0 and x._2==0
        test y._0==1 and y._1==1 and y._2==1
    end

    testset "barycentric coordinate - square" do
        terracode
            var A = hypercube.new(interval.new(1, 3), interval.new(0, 0), interval.new(0, 1))
            A:setorigin(point.from(1,0,0))
            var x = A:barycentriccoord(point.from(1,0,0))
            var y = A:barycentriccoord(point.from(3,0,1))
            var z = A:barycentriccoord(point.from(3,0,0.5))
            var p = A:barycentriccoord(point.from(1,1,0))

        end
        test x._0==0 and x._1==0 and x._2==0
        test y._0==1 and y._1==1 and y._2==0
        test z._0==1 and z._1==.5 and z._2==0
        test p._0==0 and p._1==0 and p._2==0 --projection
    end

    testset "cube evaluation - origin==center" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5), interval.new(2, 4))
        end
        test A(0,0,0) == point.from(0,1,2)
        test A(0.5,0.5,0.5) == point.from(1,3,3)
        test A(1,1,1) == point.from(2,5,4)
    end

    testset "cube evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 5), interval.new(2, 4))
            A:setorigin(point.from(2,5,4))
        end
        test A(0,0,0) == point.from(2,5,4)
        test A(1,1,1) == point.from(0,1,2)
    end

    testset "surface evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(1, 1), interval.new(2, 4))
            A:setorigin(point.from(0,1,2))
        end
        test A(0,0) == point.from(0,1,2)
        test A(1,0) == point.from(2,1,2)
    end

    testset "segment evaluation - custom origin" do
        terracode
            var A = hypercube.new(interval.new(2, 2), interval.new(1, 1), interval.new(1, 3))
            A:setorigin(point.from(2,1,3))
        end
        test A(0) == point.from(2,1,3)
        test A(1) == point.from(2,1,1)
    end

    testset "empty intersection" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(3, 4), interval.new(0, 2), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()
    end

    testset "nonempty intersection - cube" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(1, 3), interval.new(0, 2))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(0, 2), interval.new(1, 2), interval.new(0, 2))
    end

    testset "nonempty intersection - surface" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(2, 3))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(2, 2))
    end

    testset "nonempty intersection - line" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(2, 3), interval.new(0, 2), interval.new(2, 3))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(2, 2), interval.new(0, 2), interval.new(2, 2))
    end

    testset "nonempty intersection - point" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 2))
            var B = hypercube.new(interval.new(2, 3), interval.new(2, 3), interval.new(2, 3))
            var C = hypercube.intersect(A, B)
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(2, 2), interval.new(2, 2), interval.new(2, 2))
    end

    testset "product" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 0))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 0), interval.new(0, 1))
            var C = A * B
        end
        test C:isempty()==false
        test C==hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 1))
    end

    testset "division" do
        terracode
            var A = hypercube.new(interval.new(0, 2), interval.new(0, 2), interval.new(0, 0))
            var B = hypercube.new(interval.new(0, 0), interval.new(0, 0), interval.new(0, 1))
            var C = A * B
        end
        test (C / A) * A == C
        test (C / B) * B == C
    end
end
