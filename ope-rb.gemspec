$:.push File.expand_path('../lib', __FILE__)
require 'ope-rb/version'

Gem::Specification.new do |s|
  
  s.name        = 'ope-rb'
  s.version     = OPE::VERSION
  s.authors     = ['Louis Mullie']
  s.email       = ['louis.mullie@gmail.com']
  s.homepage    = 'https://github.com/cryodex/ope-rb'
  s.summary     = %q{ A Ruby implementation of order-preserving encryption }
  s.description = %q{  }

  s.files = Dir.glob('lib/**/*.rb') +
  Dir.glob('ext/**/*.{c,h,rb}')

  s.extensions << 'ext/ope/extconf.rb'
  s.add_development_dependency 'rspec', '~> 2.12.0'
  s.add_development_dependency 'rake'
  
end