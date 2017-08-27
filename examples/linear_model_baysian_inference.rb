require "rubybayes"
require "rubymc"
require "rubystats"
require "csv"


def read_csv(file)
	content = CSV.read(file)
	h = {}
	keys = content.shift		
	keys.each {|key| h[key]=[] }
	content.each do |row|
		row.each_with_index do |col,i|
			h[keys[i]] << col.to_f
		end
	end
	h
end


class ProposalDensity
  def initialize
  	@dens = []
    @dens << Rubystats::NormalDistribution.new(0.0,0.25)		
		@dens << Rubystats::NormalDistribution.new(0.0,0.01)		
		@dens << Rubystats::NormalDistribution.new(0.0,0.01)	
  end
  def rng(x)
    displacement = @dens.collect{|d| d.rng }
    x.zip(displacement).map {|x,y| x + y}
  end
end

class Posterior
  include Rubybayes::MakeLogPosterior
  include Rubybayes::EnableRandomWalk  
	
	def initialize
		data = read_csv("./howell1.csv")
		@h = data["height"]
		
		wtemp = data["weight"]
		mean = wtemp.inject(0.0) {|sum,w| sum + w} / wtemp.size
		@w = wtemp.collect{|w| w - mean}
		
		@priors = []
    @priors << Rubystats::NormalDistribution.new(178.0,100.0)
    @priors << Rubystats::NormalDistribution.new(0.0,10.0)
    @priors << Rubystats::NormalDistribution.new(1.95,0.6) 
  end
  
  def proposal_density
    ProposalDensity.new
  end

  def log_likelihood(x)
    i=0
		ll = @h.inject(0.0) do |sum,h|
			mu = x[0].to_f + x[1].to_f * @w[i]
			d = Rubystats::NormalDistribution.new(mu, Math.exp(x[2].to_f))
			i += 1
			dens = Math.log(d.pdf(h.to_f)) 			
			sum + dens
		end
    ll
  end

  def log_priors(x)
    sum = 0.0
    @priors.size.times {|i| sum + Math.log(@priors[i].pdf(x[i]))}
  end

  def starting_point
    [178.0,0.0,1.8]
  end

end


posterior = Posterior.new
puts "Model: h_i ~ N(mu,sd)"
puts "        mu ~ a + b * w_i"
puts "         a ~ N(178,100)"
puts "         b ~ N(0,10)"
puts "   log(sd) ~ N(1.95,0.67)"
puts ""
puts "        h_i = heights from Howell data"
puts "        w_i = weights from Howell data"
puts ""

# mode = posterior.mode
# puts "mode of alpha = #{mode[0]}"
# puts "mode of beta  = #{mode[1]}"
# puts "mode of sd    = #{Math.exp(mode[2])}"

10000.times { posterior.sample }
ret = []
50000.times { ret << posterior.sample }

a, b, sd = Rubymc::MonteCarloSimulation.extract_measurements({0 => ret})[0]

puts "mean of a  = #{a.mean}"
puts "mean of b  = #{b.mean}"
puts "mean of sd = #{Math.exp(sd.mean)}"
