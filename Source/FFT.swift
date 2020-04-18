//
//  FFT.swift
//  Spectra
//
//  Created by Reza Ali on 4/3/20.
//  Copyright © 2020 Hi-Rez, Inc. All rights reserved.
//

import Accelerate

open class FFT {
    public enum WindowSequence {
        case hanningNormalized
        case hanningDenormalized
        case hamming
        case blackman
        case none
    }

    public var windowSequence: FFT.WindowSequence = .blackman {
        didSet {
            setupWindow()
        }
    }

    public var smoothing: Float = 0.0

    private var samples: Int = 0
    private var n: Int = 0
    private var nHalf: Int = 0
    private var log2n: vDSP_Length = 0
    private var fft: vDSP.FFT<DSPSplitComplex>!

    private var window: [Float] = []
    private var windowedSignal: [Float] = []
    private var realParts: [Float] = []
    private var imagParts: [Float] = []
    private var spectrum: [Float] = []
    private var spectrumSmooth: [Float] = []

    public init?(samples: Int, windowSequence: FFT.WindowSequence = .blackman, smoothing: Float = 0.0) {
        let lg2 = logbf(Float(samples))
        if remainderf(Float(samples), powf(2.0, lg2)) != 0 {
            print("Number of samples is not a power of two: \(samples)")
            return nil
        }

        self.samples = samples
        self.windowSequence = windowSequence
        self.smoothing = smoothing
        setup()
    }

    private func setup() {
        n = samples
        nHalf = n / 2
        log2n = vDSP_Length(log2(Float(n)))

        guard let fft = vDSP.FFT(log2n: log2n, radix: .radix2, ofType: DSPSplitComplex.self) else {
            fatalError("Can't create FFT Setup")
        }

        self.fft = fft

        setupBuffers()
        setupWindow()
    }

    private func setupBuffers() {
        realParts = [Float](repeating: 0, count: nHalf)
        imagParts = [Float](repeating: 0, count: nHalf)
        spectrum = [Float](repeating: 0, count: nHalf)
        spectrumSmooth = [Float](repeating: 0, count: nHalf)
        windowedSignal = [Float](repeating: 0, count: n)
        window = [Float](repeating: 1, count: n)
    }

    private func setupWindow() {
        switch windowSequence {
        case .blackman:
            vDSP_blkman_window(&window, vDSP_Length(n), 0)
        case .hamming:
            vDSP_hamm_window(&window, vDSP_Length(n), 0)
        case .hanningDenormalized:
            vDSP_hann_window(&window, vDSP_Length(n), Int32(vDSP_HANN_DENORM))
        case .hanningNormalized:
            vDSP_hann_window(&window, vDSP_Length(n), Int32(vDSP_HANN_NORM))
        case .none:
            window = [Float](repeating: 1, count: n)
        }
    }

    public func forward(_ signal: UnsafeMutablePointer<Float>) {
        // Window Signal
        window.withUnsafeMutableBufferPointer { windowPtr in
            windowedSignal.withUnsafeMutableBufferPointer { windowedSignalPtr in
                vDSP_vmul(signal, 1, windowPtr.baseAddress!, 1, windowedSignalPtr.baseAddress!, 1, vDSP_Length(n))
            }
        }

        _forward()
    }

    public func forward(_ signal: [Float]) {
        // Window Signal
        window.withUnsafeMutableBufferPointer { windowPtr in
            windowedSignal.withUnsafeMutableBufferPointer { windowedSignalPtr in
                vDSP_vmul(signal, 1, windowPtr.baseAddress!, 1, windowedSignalPtr.baseAddress!, 1, vDSP_Length(n))
            }
        }

        _forward()
    }

    private func _forward() {
        // Perform FFT
        spectrum.withUnsafeMutableBufferPointer { specPtr in
            realParts.withUnsafeMutableBufferPointer { realPtr in
                imagParts.withUnsafeMutableBufferPointer { imagPtr in

                    // Create a `DSPSplitComplex` to contain the signal.
                    var complexSignal = DSPSplitComplex(realp: realPtr.baseAddress!,
                                                        imagp: imagPtr.baseAddress!)

                    // Convert to complex numbers
                    windowedSignal.withUnsafeBytes {
                        vDSP.convert(interleavedComplexVector: [DSPComplex]($0.bindMemory(to: DSPComplex.self)),
                                     toSplitComplexVector: &complexSignal)
                    }

                    let length = vDSP_Length(nHalf)
                    // Perform FFT
                    fft.forward(input: complexSignal, output: &complexSignal)
                    
                    // Process signal square root of the absolute value of each element of imagParts.
                    vDSP_zvmags(&complexSignal, 1, specPtr.baseAddress!, 1, length)

                    // vDSP_zvmagsD returns squares of the FFT magnitudes, so take the root here
                    var len = Int32(nHalf)
                    vvsqrtf(specPtr.baseAddress!, specPtr.baseAddress!, &len)

                    // Scale to get into a good 0 - 1 range
                    var scalar: [Float] = [8.0 / Float(nHalf)]
                    vDSP_vsmul(specPtr.baseAddress!, 1, &scalar, specPtr.baseAddress!, 1, vDSP_Length(nHalf))
                }
            }
        }

        if smoothing > 0.0 {
            spectrum.withUnsafeMutableBufferPointer { specPtr in
                spectrumSmooth.withUnsafeMutableBufferPointer { specSmoothPtr in
                    vDSP_vsmul(specSmoothPtr.baseAddress!, 1, &smoothing, specSmoothPtr.baseAddress!, 1, vDSP_Length(nHalf))
                    vDSP_vmaxmg(specPtr.baseAddress!, 1, specSmoothPtr.baseAddress!, 1, specSmoothPtr.baseAddress!, 1, vDSP_Length(nHalf))
                }
            }
        }
    }

    public func getSpectrum() -> [Float] {
        if smoothing > 0.0 {
            return spectrumSmooth
        }
        return spectrum
    }

    public func getReal() -> [Float] {
        return realParts
    }

    public func getImaginary() -> [Float] {
        return imagParts
    }

    public func getWindow() -> [Float] {
        return window
    }

    deinit {        
    }
}
