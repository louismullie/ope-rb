$:.push File.expand_path('../lib', __FILE__)
require 'ope/version'

Gem::Specification.new do |s|
  
  s.name        = 'ope'
  s.version     = OPE::VERSION
  s.authors     = ['Louis Mullie']
  s.email       = ['louis.mullie@gmail.com']
  s.homepage    = 'https://github.com/symeapp/srp'
  s.summary     = %q{ Fast C implementation of order-preserving encryption. }
  s.description = %q{  }

  s.files = Dir.glob('lib/**/*.rb') +
  Dir.glob('ext/**/*.{c,h,rb}')

  s.extensions << 'ext/ope/extconf.rb'
  s.add_development_dependency 'rspec', '~> 2.12.0'
  s.add_development_dependency 'rake'
  
end