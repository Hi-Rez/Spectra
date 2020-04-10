//
//  FFT.swift
//  Spectra
//
//  Created by Reza Ali on 4/3/20.
//  Copyright © 2020 Hi-Rez, Inc. All rights reserved.
//

import Accelerate

open class FFT {
    public var windowSequence: vDSP.WindowSequence = .blackman {
        didSet {
            setupWindow()
        }
    }

    public var useWindow: Bool = true {
        didSet {
            if !useWindow {
                window = [Float](repeating: 1, count: n)
            }
        }
    }

    public var samples: Int = 0 {
        didSet {
            setup()
        }
    }

    private var n: Int = 0
    private var nHalf: Int = 0
    private var log2n: vDSP_Length = 0
    private var fft: vDSP.FFT<DSPSplitComplex>!

    private var window: [Float] = []
    private var windowedSignal: [Float] = []
    private var realParts: [Float] = []
    private var imagParts: [Float] = []
    private var spectrum: [Float] = []

    public init(samples: Int, windowSequence: vDSP.WindowSequence = .blackman) {
        self.samples = samples
        self.windowSequence = windowSequence
        setup()
    }

    private func setup() {
        let lg2 = logbf(Float(samples))
        if remainderf(Float(samples), powf(2.0, lg2)) != 0 {
            print("Number of samples is not a power of two: \(samples)")
            return
        }

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
        windowedSignal = [Float](repeating: 0, count: n)
        window = [Float](repeating: 1, count: n)
    }

    private func setupWindow() {
        if useWindow {
            window = vDSP.window(ofType: Float.self, usingSequence: windowSequence, count: n, isHalfWindow: false)
        }
        else
        {
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

                    // Perform FFT
                    fft.forward(input: complexSignal, output: &complexSignal)

                    // Process signal square root of the absolute value of each element of imagParts.
                    vDSP_zvmags(&complexSignal, 1, specPtr.baseAddress!, 1, vDSP_Length(nHalf))

                    // vDSP_zvmagsD returns squares of the FFT magnitudes, so take the root here
                    var len = Int32(nHalf)
                    vvsqrtf(specPtr.baseAddress!, specPtr.baseAddress!, &len)

                    // Scale to get into a good 0 - 1 range
                    var scalar: [Float] = [2.0 / Float(nHalf)]
                    vDSP_vsmul(specPtr.baseAddress!, 1, &scalar, specPtr.baseAddress!, 1, vDSP_Length(nHalf))
                }
            }
        }
    }


    public func forward(_ signal: [Float]) {
        // Window Signal
        window.withUnsafeMutableBufferPointer { windowPtr in
            windowedSignal.withUnsafeMutableBufferPointer { windowedSignalPtr in
                vDSP_vmul(signal, 1, windowPtr.baseAddress!, 1, windowedSignalPtr.baseAddress!, 1, vDSP_Length(n))
            }
        }

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

                    // Perform FFT
                    fft.forward(input: complexSignal, output: &complexSignal)

                    // Process signal square root of the absolute value of each element of imagParts.
                    vDSP_zvmags(&complexSignal, 1, specPtr.baseAddress!, 1, vDSP_Length(nHalf))

                    // vDSP_zvmagsD returns squares of the FFT magnitudes, so take the root here
                    var len = Int32(nHalf)
                    vvsqrtf(specPtr.baseAddress!, specPtr.baseAddress!, &len)

                    // Scale to get into a good 0 - 1 range
                    var scalar: [Float] = [2.0 / Float(nHalf)]
                    vDSP_vsmul(specPtr.baseAddress!, 1, &scalar, specPtr.baseAddress!, 1, vDSP_Length(nHalf))
                }
            }
        }
    }

    public func getSpectrum() -> [Float] {
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
        print("Destroying FFT")
    }
}
