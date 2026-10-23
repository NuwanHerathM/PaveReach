@testset "Volume" begin
    @testset "Single box / interval" begin
        segment = interval(-2.0, 4.5)
        box = IntervalBox(interval(-1.0, 2.0), interval(-40.0, -38.0))
        @testset "Standard" begin
            @test volume_box(segment) == 6.5
            @test volume_box(box) == 6.0
        end
        @testset "With scaling" begin
            @test volume_box(segment, 2) == 13.0
            @test volume_box(segment, [2]) == 13.0
            @test_throws DimensionMismatch volume_box(segment, [2, 3])
            @test_throws DimensionMismatch volume_box(box, 2)
            @test_throws DimensionMismatch volume_box(box, [2])
            @test volume_box(box, 0.5, 1.5) == 4.5
            @test volume_box(box, [0.5, 1.5]) == 4.5
        end
    end
    @testset "List of boxes / intervals" begin
        segments = [interval(-2.0, 4.5), interval(10.0, 11.5)]
        boxes = [IntervalBox(interval(-1.0, 2.0), interval(-40.0, -38.0)), IntervalBox(interval(0.0, 3.0), interval(-20.0, -10.0))]
        @testset "Standard" begin
            @test volume_boxes(segments) == 8.0
            @test volume_boxes(boxes) == 36.0
        end
        @testset "With scaling" begin
            @test volume_boxes(segments, 2) == 16.0
            @test volume_boxes(segments, [2]) == 16.0
            @test_throws DimensionMismatch volume_boxes(segments, [2, 3])
            @test_throws DimensionMismatch volume_boxes(boxes, 2)
            @test_throws DimensionMismatch volume_boxes(boxes, [2])
            @test volume_boxes(boxes, 0.5, 1.5) == 27.0
            @test volume_boxes(boxes, [0.5, 1.5]) == 27.0
        end
    end
end