using Base.Test
using RoboticsToolbox.transform

@testset "Rotation Tests" begin
    @testset "2D" begin
        @test ≈(rot2(0.3), [0.9553 -0.2955;
                            0.2955 0.9553],
                atol=1e-4)
        @test ≈(se2(1, 2, 0.3), [0.9553 -0.2955 1;
                                 0.2955 0.9553 2;
                                 0 0 1],
                atol=1e-4)
        @test ≈(trot2(0.3), [0.9553 -0.2955 0;
                             0.2955 0.9553 0;
                             0 0 1],
                atol=1e-4)
    end

    @testset "3D" begin
        @testset "Rotation Matrices" begin
            @testset "No rotation" begin
                @test ≈(rotx(0), eye(3), atol=1e-9)
                @test ≈(roty(0), eye(3), atol=1e-9)
                @test ≈(rotz(0), eye(3), atol=1e-9)
            end

            @testset "Rotation by π/2" begin
                @test ≈(rotx(π / 2), [1 0 0;
                                      0 0 -1;
                                      0 1 0],
                        atol=1e-9)
                @test ≈(roty(π / 2), [0 0 1;
                                      0 1 0;
                                      -1 0 0],
                        atol=1e-9)
                @test ≈(rotz(π / 2), [0 -1 0;
                                      1 0 0;
                                      0 0 1],
                        atol=1e-9)
            end
        end

        @testset "Homogeneous Transforms" begin
            @testset "No rotation" begin
                @test ≈(trotx(0), eye(4), atol=1e-9)
                @test ≈(troty(0), eye(4), atol=1e-9)
                @test ≈(trotz(0), eye(4), atol=1e-9)
            end

            @testset "Rotation by π/2" begin
                @test ≈(trotx(π / 2), [1 0 0 0;
                                       0 0 -1 0;
                                       0 1 0 0;
                                       0 0 0 1],
                        atol=1e-9)
                @test ≈(troty(π / 2), [0 0 1 0;
                                       0 1 0 0;
                                       -1 0 0 0;
                                       0 0 0 1],
                        atol=1e-9)
                @test ≈(trotz(π / 2), [0 -1 0 0;
                                       1 0 0 0;
                                       0 0 1 0;
                                       0 0 0 1],
                        atol=1e-9)
            end

            @testset "Rotation by 0.1" begin
                @test ≈(trotx(0.1), [1 0 0 0;
                                     0 0.995 -0.0998 0;
                                     0 0.0998 0.995 0;
                                     0 0 0 1],
                        atol=1e-4)
                @test ≈(troty(0.1), [0.995 0 0.0998 0;
                                     0 1 0 0;
                                     -0.0998 0 0.995 0;
                                    0 0 0 1],
                        atol=1e-4)
                @test ≈(trotz(0.1), [0.995 -0.0998 0 0;
                                     0.0998 0.995 0 0;
                                     0 0 1 0;
                                     0 0 0 1],
                        atol=1e-4)
            end
        end

        @test ≈(se3(se2(1, 2, 0.3)), [0.9553 -0.2955 0 1;
                                      0.2955 0.9553 0 2;
                                      0 0 1 0;
                                      0 0 0 1],
                atol=1e-4)
    end
end

@testset "Conversions" begin
    @testset "r2t" begin
        @test ≈(r2t([1 2; 3 4]), [1 2 0;
                                  3 4 0;
                                  0 0 1],
                atol=1e-4)
        @test ≈(r2t([1 2 3;
                     4 5 6;
                     7 8 9]),
                    [1 2 3 0;
                     4 5 6 0;
                     7 8 9 0;
                     0 0 0 1],
                    atol=1e-4)
    end

    @testset "t2r" begin
        @test ≈(t2r([1 2 0;
                     3 4 0;
                     0 0 1]),
                [1 2; 3 4],
                atol=1e-4)
        @test ≈(t2r([1 2 3 0;
                     4 5 6 0;
                     7 8 9 0;
                     0 0 0 1]),
                [1 2 3;
                 4 5 6;
                 7 8 9],
                atol=1e-4)
    end

    @testset "Jacobian" begin
        @test ≈(rpy2jac(0.1, 0.2, 0.3), [1 0 0.1987;
                                         0 0.995 -0.0978;
                                         0 0.0998 0.9752],
                atol=1e-4)
        @test ≈(rpy2jac(0, 0, 0), eye(3), atol=1e-4)
    end

    @testset "Roll-pitch-yaw" begin
        @testset "rpy2r" begin
            @test ≈(rpy2r(0.1, 0.2, 0.3), [0.9363 -0.2896 0.1987;
                                           0.3130 0.9447 -0.0978;
                                           -0.1593 0.1538 0.9752],
                    atol=1e-4)
            @test ≈(rpy2r(0, 0, 0), eye(3), atol=1e-4)
        end

        @testset "rpy2t" begin
            @test ≈(rpy2t(0.1, 0.2, 0.3), [0.9363 -0.2896 0.1987 0;
                                           0.3130 0.9447 -0.0978 0;
                                           -0.1593 0.1538 0.9752 0;
                                           0 0 0 1],
                    atol=1e-4)
            @test ≈(rpy2t(0, 0, 0), [eye(3) [0; 0; 0];
                                     0 0 0 1],
                    atol=1e-4)
        end

        @testset "tr2rpy" begin
            @test ≈(collect(tr2rpy(rpy2t(0.1, 0.2, 0.3))),
                    collect((0.1, 0.2, 0.3)), atol=1e-4)
        end
    end
end
