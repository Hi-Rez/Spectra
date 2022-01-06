// swift-tools-version:5.5

import PackageDescription

let package = Package(
    name: "Spectra",
    platforms: [.macOS(.v10_15), .iOS(.v13), .tvOS(.v13)],
    products: [
        .library(name: "Spectra", targets: ["Spectra"])
    ],
    targets: [
        .target(name: "Spectra")
    ],
    swiftLanguageVersions: [.v5]
)
