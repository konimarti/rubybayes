# coding: utf-8
lib = File.expand_path("../lib", __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require "rubybayes/version"

Gem::Specification.new do |spec|
  spec.name          = "rubybayes"
  spec.version       = Rubybayes::VERSION
  spec.authors       = ["konimarti"]
  spec.email         = ["koni.marti@gmail.com"]

  spec.summary       = %q{Rubybayes provides the functionality to perform complex Baysian inference by means of Monte Carlo sampling from a posterior density.}
  spec.description   = %q{Rubybayes comes with a Metropolis-Hastings Monte Carlo sampler and analysis tools to evaluate Baysian posterior densities. }
  spec.homepage      = ""
  spec.license       = "MIT"

  # Prevent pushing this gem to RubyGems.org. To allow pushes either set the 'allowed_push_host'
  # to allow pushing to a single host or delete this section to allow pushing to any host.
  if spec.respond_to?(:metadata)
    spec.metadata["allowed_push_host"] = "TODO: Set to 'http://mygemserver.com'"
  else
    raise "RubyGems 2.0 or newer is required to protect against " \
      "public gem pushes."
  end

  spec.files         = `git ls-files -z`.split("\x0").reject do |f|
    f.match(%r{^(test|spec|features)/})
  end
  spec.bindir        = "exe"
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.15"
  spec.add_development_dependency "rake", "~> 10.0"
  spec.add_development_dependency "minitest", "~> 5.0"
  spec.add_development_dependency "rubystats", "~> 0.2.6"
  spec.add_development_dependency "minimization", "~> 0.2.1"
end

