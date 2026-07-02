include("../src/pave.jl")

box = IntervalBox(interval(0, 2), interval(40, 50))
@testset "Bisect increasing / decreasing order" begin
    @testset "Dim 1" begin
        dim = 1
        box_1, box_2 = bisect(box, dim)
        @test bisect_increasing_order(box, dim) == (box_1, box_2)
        @test bisect_decreasing_order(box, dim) == (box_2, box_1)
    end
    @testset "Dim 2" begin
        dim = 2
        box_1, box_2 = bisect(box, dim)
        @test bisect_increasing_order(box, dim) == (box_1, box_2)
        @test bisect_decreasing_order(box, dim) == (box_2, box_1)
    end
end

@testset "Create bisect in optimization direction" begin
    @testset "Positive sign" begin
        sign = 1
        bisect_in_optimization_direction = create_bisect_in_optimization_direction(sign)
        dim = 1
        box_1, box_2 = bisect(box, dim)
        @test bisect_in_optimization_direction(box, dim) == (box_1, box_2)
    end
    @testset "Negative sign" begin
        sign = -1
        bisect_in_optimization_direction = create_bisect_in_optimization_direction(sign)
        dim = 1
        box_1, box_2 = bisect(box, dim)
        @test bisect_in_optimization_direction(box, dim) == (box_2, box_1)
    end
end

@testset "Bisect on precision" begin
    ϵ = [1.0, 15]
    @test IntervalArithmetic.diam(box[1]) > ϵ[1]
    @test IntervalArithmetic.diam(box[2]) < ϵ[2]
    box_1, box_2 = bisect(box, 1)
    @test bisect_precision(box, ϵ) == (box_1, box_2, 1)
end

@testset "Slices according to optimization" begin
    a = interval(0, 5)
    @test backward_bound(a, -1) == 5
    @test forward_bound(a, -1) == 0
    @test backward_bound(a, 1) == 0
    @test forward_bound(a, 1) == 5
    optimization_directions = [-1, 1]
    @test backward_slices(box, optimization_directions, nothing) == [IntervalBox(interval(2, 2), interval(40, 50)), IntervalBox(interval(0, 2), interval(40, 40))]
    @test forward_slices(box, optimization_directions, nothing) == [IntervalBox(interval(0, 0), interval(40, 50)), IntervalBox(interval(0, 2), interval(50, 50))]
    unchanged_face = (2, Lower)
    @test backward_slices(box, optimization_directions, unchanged_face) == [IntervalBox(interval(2, 2), interval(40, 50))]
    @test forward_slices(box, optimization_directions, unchanged_face) == [IntervalBox(interval(0, 0), interval(40, 50)), IntervalBox(interval(0, 2), interval(50, 50))]
end

@testset "Atom conjunction" begin
    atom = Atom("P_1")
    @test atom_conjunction(atom, 1) == atom
    @test atom_conjunction(atom, 2) == parseformula("P_1_1 ∧ P_1_2")
    @test atom_conjunction(atom, 4) == parseformula("((P_1_1 ∧ P_1_2) ∧ P_1_3) ∧ P_1_4")
    @test_throws MethodError atom_conjunction(:SyntaxBranch, :Int)
end

@testset "Atom complement disjunction" begin
    atom = Atom("P")
    intervals = [interval(0,1)]
    @test atom_complement_disjunction(atom, intervals) == parseformula("P_1 ∨ P_2")
    intervals = [interval(minus_inf, 4)]
    @test atom_complement_disjunction(atom, intervals) == parseformula("P")
    intervals = [interval(0, 1), interval(minus_inf, 4), interval(7, 8)]
    @test atom_complement_disjunction(atom, intervals) == parseformula("((P_1_1 ∨ P_1_2) ∨ P_2) ∨ (P_3_1 ∨ P_3_2)")
    @test_throws MethodError atom_complement_disjunction(:SyntaxBranch, :Int)
end

@testset "Expand" begin
    M = parseformula("(¬ P_1 ∧ P_2) ∨ P_3")
    intervals = [interval(-5, 5), interval(-5, 5), interval(minus_inf, 8), interval(3, plus_inf), interval(0, 10), interval(1, 2)]
    sizes = [2, 3, 1]
    indices_dict = build_ranges_dict(M, sizes)
    @test expand_1(M, indices_dict, intervals) == parseformula("(((P_1_1_1 ∨ P_1_1_2) ∨ (P_1_2_1 ∨ P_1_2_2)) ∧ (P_2_1 ∧ P_2_2 ∧ P_2_3)) ∨ P_3")
    @test expand_2(M, indices_dict, intervals) == parseformula("((P_1_1 ∧ P_1_2) ∨ (P_2_1 ∨ P_2_2 ∨ (P_2_3_1 ∨ P_2_3_2))) ∧ (P_3_1 ∨ P_3_2)")
    M = parseformula("P_1 ∧ ¬ (P_2 ∨ P_3)")
    @test_throws ErrorException expand_1(M, indices_dict, intervals)
    @test_throws ErrorException expand_2(M, indices_dict, intervals)
end

# @testset "Expand negations" begin
#     M = parseformula("P_1")
#     intervals = [interval(0, 6)]
#     sizes = [1]
#     indices_dict = build_ranges_dict(M, sizes)
#     @test expand_negations_in_1(M, indices_dict, intervals) == parseformula("P_1")
#     @test expand_negations_out_1(M, indices_dict, intervals) == parseformula("P_1")
#     @test expand_negations_in_2(M, indices_dict, intervals) == parseformula("P_1_1 ∨ P_1_2")
#     @test expand_negations_out_2(M, indices_dict, intervals) == parseformula("P_1_1 ∨ P_1_2")
#     M = parseformula("¬ P_1")
#     intervals = [interval(0, 6)]
#     sizes = [1]la
#     indices_dict = build_ranges_dict(M, sizes)
#     @test expand_negations_in_1(M, indices_dict, intervals) == parseformula("P_1_1 ∨ P_1_2")
#     @test expand_negations_out_1(M, indices_dict, intervals) == parseformula("P_1_1 ∨ P_1_2")
#     @test expand_negations_in_2(M, indices_dict, intervals) == parseformula("P_1")
#     @test expand_negations_out_2(M, indices_dict, intervals) == parseformula("P_1")
#     M = parseformula("P_1 ∨ (¬ P_2 ∧ ¬ P_3) ∧ (P_4 ∨ P_5) ∧ ¬ P_6")
#     intervals = [interval(-5, 5), interval(-5, 5), interval(minus_inf, 8), interval(3, plus_inf), interval(0, 10), interval(1, 2)]
#     sizes = repeat([1], 6)
#     indices_dict = build_ranges_dict(M, sizes)
#     @test expand_negations_in_1(M, indices_dict, intervals) == parseformula("P_1 ∨ ((P_2_1 ∨ P_2_2) ∧ P_3) ∧ (P_4 ∨ P_5) ∧ (P_6_1 ∨ P_6_2)")
#     @test expand_negations_out_1(M, indices_dict, intervals) == parseformula("P_1 ∨ ((P_2_1 ∨ P_2_2) ∧ P_3) ∧ (P_4 ∨ P_5) ∧ (P_6_1 ∨ P_6_2)")
#     @test expand_negations_in_2(M, indices_dict, intervals) == parseformula("(P_1_1 ∨ P_1_2) ∧ ((P_2 ∨ P_3) ∨ (P_4 ∧ (P_5_1 ∨ P_5_2)) ∨ P_6)")
#     @test expand_negations_out_2(M, indices_dict, intervals) == parseformula("(P_1_1 ∨ P_1_2) ∧ ((P_2 ∨ P_3) ∨ (P_4 ∧ (P_5_1 ∨ P_5_2)) ∨ P_6)")
# end

@testset "Positions / ranges" begin
    @testset "Scalar functions" begin
        M = parseformula("P_1 ∨ (¬ P_2 ∧ ¬ P_3) ∧ (P_4 ∨ P_5) ∧ ¬ P_6")
        intervals = [interval(-5, 5), interval(-5, 5), interval(minus_inf, 8), interval(3, plus_inf), interval(0, 10), interval(1, 2)]
        sizes = repeat([1], 6)
        indices_dict = build_ranges_dict(M, sizes)
        @testset "Duplication positions" begin
            @test issetequal(duplication_positions_1(M, indices_dict, intervals), [2, 6])
            @test issetequal(duplication_positions_2(M, indices_dict, intervals), [1, 5])
        end
        @testset "Complement ranges" begin
            @test issetequal(complement_ranges_1(M, indices_dict), [2:2, 3:3, 6:6])
            @test issetequal(complement_ranges_2(M, indices_dict), [1:1, 4:4, 5:5])
        end
    end
    @testset "Vectorial functions" begin
        M = parseformula("(¬ P_1 ∧ P_2) ∨ P_3")
        intervals = [interval(-5, 5), interval(-5, 5), interval(minus_inf, 8), interval(3, plus_inf), interval(0, 10), interval(1, 2)]
        sizes = [2, 3, 1]
        indices_dict = build_ranges_dict(M, sizes)
        @testset "Duplication positions" begin
            @test issetequal(duplication_positions_1(M, indices_dict, intervals), [1, 2])
            @test issetequal(duplication_positions_2(M, indices_dict, intervals), [5, 6])
        end
        @testset "Complement ranges" begin
            @test issetequal(complement_ranges_1(M, indices_dict), [1:2])
            @test issetequal(complement_ranges_2(M, indices_dict), [3:5, 6:6])
        end
    end
end

@testset "Expand intervals" begin
    intervals = [interval(-5, 5), interval(-5, 5), interval(minus_inf, 8), interval(3, plus_inf), interval(0, 10), interval(1, 2)]
    # sizes = repeat([1], 6)
    @test expand_intervals(intervals, [2:2, 6:6]) == [interval(-5, 5), interval(minus_inf, -5 - strict_epsilon), interval(5 + strict_epsilon, plus_inf), interval(minus_inf, 8), interval(3, plus_inf), interval(0, 10), interval(minus_inf, 1 - strict_epsilon), interval(2 + strict_epsilon, plus_inf)]
    @test expand_intervals(intervals, [1:1, 5:5]) == [interval(minus_inf, -5 - strict_epsilon), interval(5 + strict_epsilon, plus_inf), interval(-5, 5), interval(minus_inf, 8), interval(3, plus_inf), interval(minus_inf, -strict_epsilon), interval(10 + strict_epsilon, plus_inf), interval(1, 2)]
    # sizes = [2, 3, 1]
    @test expand_intervals(intervals, [1:2, 6:6]) == [interval(minus_inf, -5 - strict_epsilon), interval(5 + strict_epsilon, plus_inf), interval(minus_inf, -5 - strict_epsilon), interval(5 + strict_epsilon, plus_inf), interval(minus_inf, 8), interval(3, plus_inf), interval(0, 10), interval(minus_inf, 1 - strict_epsilon), interval(2 + strict_epsilon, plus_inf)]
end

@testset "Inflate list" begin
    @test inflate([7,6,3,5,0], [2,4]) == [7,6,6,3,5,5,0]
end

@testset "Zero membership" begin
    R = [EmptySet(1), LazySets.IntervalModule.Interval(4, 6), LazySets.IntervalModule.Interval(-100, 5)]
    @test is_zero_in(R, 1) == false
    @test is_zero_in(R, 2) == false
    @test is_zero_in(R, 3) == true
    @test is_zero_not_in(R, 1) == true
    @test is_zero_not_in(R, 2) == true
    @test is_zero_not_in(R, 3) == false
end

@testset "Oracles" begin
    M = parseformula("P_1 ∧ ((P_2 ∨ P_3) ∧ P_4)")
    R = [EmptySet(1), LazySets.IntervalModule.Interval(4, 6), LazySets.IntervalModule.Interval(-100, 5), LazySets.IntervalModule.Interval(-2, 0.5)]
    positions_dict = Dict(Atom("P_1") => 1, Atom("P_2") => 2, Atom("P_3") => 3, Atom("P_4") => 4)
    @test test_zero_in(M, R, positions_dict) == false
    @test test_zero_not_in(M, R, positions_dict) == true
    R = [LazySets.IntervalModule.Interval(-10, 10), LazySets.IntervalModule.Interval(4, 6), LazySets.IntervalModule.Interval(-100, 5), LazySets.IntervalModule.Interval(-2, 0.5)]
    @test test_zero_in(M, R, positions_dict) == true
    @test test_zero_not_in(M, R, positions_dict) == false
end