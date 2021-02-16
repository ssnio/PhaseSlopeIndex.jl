using PhaseSlopeIndex
using Test
using Statistics: mean

@testset "PhaseSlopeIndex.jl" begin
    # tests of data2psi ###############################################
    # two random signals
    signal = [[rand(-1.0:0.0001:1.0, 1000000);] [rand(-1.0:0.0001:1.0, 1000000);]]
    psi = data2psi(signal, 100)[1]
    @test all(isapprox(psi, zeros(2, 2); atol=0.01))

    # channel 2 being exactly channel 1
    ch1_ = rand(100000)
    signal = [[ch1_;] [ch1_;]]
    psi = data2psi(signal, 100)[1]
    @test psi == zeros(2, 2)

    # channel 1 leading channel 2
    ch1_ = rand(100000)
    signal = [[ch1_[2:end];] [ch1_[1:(end - 1)];]]
    psi = data2psi(signal, 100)[1]
    @test psi[1, 1] == 0.0
    @test psi[2, 2] == 0.0
    @test psi[1, 2] > 1.0
    @test psi[2, 1] + psi[1, 2] == 0

    # channel 2 leading channel 1
    ch2_ = rand(100000)
    signal = [[ch2_[1:(end - 1)];] [ch2_[2:end];]]
    psi = data2psi(signal, 100)[1]
    @test psi[1, 1] == 0.0
    @test psi[2, 2] == 0.0
    @test psi[1, 2] < -1.0
    @test psi[2, 1] + psi[1, 2] == 0

    # test data2para ##################################################
    signal = rand(100, 3, 2)
    @test_throws DimensionMismatch data2psi(signal, 100)

    # tests of int ####################################################
    @test PhaseSlopeIndex.int(3.14) == 3
    @test PhaseSlopeIndex.int(-2.72) == -2
    @test PhaseSlopeIndex.int(0.0) == 0

    # tests of dropmean ###############################################
    test_a = rand(13)
    @test PhaseSlopeIndex.dropmean(test_a, 1) == mean(test_a; dims=1)
    test_a = rand(3, 5)
    @test PhaseSlopeIndex.dropmean(test_a, 1) == mean(test_a; dims=1)[1, :]
    test_a = rand(7, 11)
    @test PhaseSlopeIndex.dropmean(test_a, 2) == mean(test_a; dims=2)[:, 1]

    # tests of squeeze ################################################
    test_a = rand(17)
    @test PhaseSlopeIndex.squeeze(test_a) == test_a
    test_a = rand(3, 5, 7)
    @test PhaseSlopeIndex.squeeze(test_a) == test_a
    test_a = rand(3, 1, 7)
    @test PhaseSlopeIndex.squeeze(test_a) == test_a[:, 1, :]

    # tests of detrend ################################################
    # testing the inplace
    test_x = [0:0.01:2.0;]
    zeros_x = zeros(size(test_x))
    @test all(isapprox.(PhaseSlopeIndex.detrend!(test_x, 1), zeros_x, atol=0.001))
    @test all(isapprox.(test_x, zeros_x, atol=0.001))

    test_x = [0:0.01:2.0;]
    @test all(isapprox.(PhaseSlopeIndex.detrend!(test_x, 0), [-1.0:0.01:1.0;], atol=0.001))
    test_x = [0:0.01:2.0;]
    @test all(isapprox.(PhaseSlopeIndex.detrend!(test_x, 1), zeros_x, atol=0.001))

    test_x = [0:0.01:(2pi);]
    test_y = sin.(test_x)
    @test all(isapprox.(PhaseSlopeIndex.detrend!(test_y, 0), test_y, atol=0.001))
    test_x = [0:0.01:(2pi);]
    test_y = sin.(test_x)
    @test all(isapprox.(PhaseSlopeIndex.detrend!(test_y, 1), test_y, atol=0.001))

    test_x = [0:0.01:(2pi);]
    test_y = cos.(test_x) .+ test_x
    aux_x = [0:0.01:(2pi);] .- pi
    aux_y = cos.(test_x) .+ aux_x
    @test all(isapprox.(PhaseSlopeIndex.detrend!(test_y, 0), aux_y, atol=0.001))
    test_x = [0:0.01:(2pi);]
    test_y = cos.(test_x) .+ test_x
    @test all(isapprox.(PhaseSlopeIndex.detrend!(test_y, 1), cos.(test_x), atol=0.01))

    # tests of hanning_fun ############################################
    # the arrays are from MATLAB `hanning` function.
    hann_12 = [
        0.057271987173395
        0.215967626634422
        0.439731659872338
        0.677302443521268
        0.874255374085551
        0.985470908713026
        0.985470908713026
        0.874255374085551
        0.677302443521268
        0.439731659872338
        0.215967626634422
        0.057271987173395
    ]
    hann_13 = [
        0.049515566048790
        0.188255099070633
        0.388739533021843
        0.611260466978157
        0.811744900929367
        0.950484433951210
        1.000000000000000
        0.950484433951210
        0.811744900929367
        0.611260466978157
        0.388739533021843
        0.188255099070633
        0.049515566048790
    ]
    hanning_window_12 = PhaseSlopeIndex.hanning_fun(12)
    hanning_window_13 = PhaseSlopeIndex.hanning_fun(13)
    @test all(hanning_window_12 .≈ hann_12)
    @test all(hanning_window_13 .≈ hann_13)
    @test hanning_window_12[1] == hanning_window_12[end]
    @test hanning_window_13[1] == hanning_window_13[end]
    @test hanning_window_12[1] > 0.0
    @test hanning_window_13[1] > 0.0
end
