Pod::Spec.new do |spec|
  spec.name                   = "Spectra"
  spec.version                = "0.0.1"
  spec.summary                = "Spectra helps to analyze a signal and get it spectrum via Apple's Accelerate FFT functions"
  spec.description            = <<-DESC
  Spectra helps to analyze a signal and get it FFT spectrum via Apple's Accelerate framework
                   DESC
  spec.homepage               = "https://github.com/Hi-Rez/Spectra"
  spec.license                = { :type => "MIT", :file => "LICENSE" }
  spec.author                 = { "Reza Ali" => "reza@hi-rez.io" }
  spec.social_media_url       = "https://twitter.com/rezaali"
  spec.source                 = { :git => "https://github.com/Hi-Rez/Spectra.git", :tag => spec.version.to_s }

  spec.osx.deployment_target  = "10.10"
  spec.ios.deployment_target  = "4.0"
  spec.tvos.deployment_target = "9.0"

  spec.source_files           = "Source/*.h", "Source/**/*.{h,m,swift}"
  spec.frameworks             = "Accelerate"
  spec.module_name            = "Spectra"
  spec.swift_version          = "5.1"
end

