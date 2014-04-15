require 'rspec/core/rake_task'
require 'rake/clean'
require 'rspec/core/rake_task'

desc 'run specs'
task :spec do
  RSpec::Core::RakeTask.new
end

desc 'clean, compile C ext and copy objects into lib'
task :make => [:clean] do
  Dir.chdir('ext/ope') do
    ruby 'extconf.rb'
    sh 'make'
  end
  cp 'ext/ope/native.bundle', 'lib/ope'
end

CLEAN.include('ext/**/*{.o,.log,.so,.bundle}')
CLEAN.include('lib/**/*{.o,.log,.so,.bundle}')
CLEAN.include('ext/**/Makefile')