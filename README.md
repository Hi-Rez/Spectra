# Spectra

## About

Spectra helps to analyze a signal and get it FFT spectrum via Apple's Accelerate framework

## Example

```
let numberOfSamples: Int = 512 // must be a power of two

//Do this once when setting up or when your number of samples change
var fft = FFT(samples: numberOfSamples)

var buffer: [Float] = [...]
fft.forward(buffer)

var spectrum: [Float] = fft.getSpectrum()
```
