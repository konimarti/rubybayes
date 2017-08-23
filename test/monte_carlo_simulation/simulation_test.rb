require "test_helper"
require "rubystats"

class SimulationTest < Minitest::Test
  def test_simulation    
    experiment = Rubybayes::MonteCarloSimulation::Simulation.new do
      burn_in 10
      iterations 100
      sample { 1.0 }
      calculate {|x| x * 2.0 }
    end    
    
    results = experiment.run[0]      
    assert_equal results.size, 100
    
    mean = results.inject(0.0) {|sum,x| sum + x} / results.size
    assert_in_epsilon 2.0, mean, 0.000001   
  end
  
  def test_simulation_short
    results = Rubybayes::MonteCarloSimulation.run do
      burn_in 10
      iterations 100
      sample { 1.0 }
      calculate {|x| x * 2.0 }
    end    
      
    assert_equal results[0].size, 100
    
    mean = results[0].inject(0.0) {|sum,x| sum + x} / results[0].size
    assert_in_epsilon 2.0, mean, 0.000001   
  end
  
  def test_calculate_pi_with_multiple_chains
    results = Rubybayes::MonteCarloSimulation.run do
      burn_in 0
      iterations 4000
      chains 5
      sample { [Kernel.rand, Kernel.rand] }
      calculate {|x| (Math.sqrt(x[0]**2+x[1]**2) <= 1.0) ? 1.0 : 0.0 }
    end    
    
    results.each_pair do |chain, result|
      assert_equal result.size, 4000
    end   
    
    combined = Rubybayes::MonteCarloSimulation.merge_chains(results)
    
    m = Rubybayes::MonteCarloSimulation::Measurement.new(combined[0])
    
    assert_in_epsilon 4.0*m.mean, Math::PI, 0.1
    
  end
  
  def test_simulation_short_multiple_dimensions
    results = Rubybayes::MonteCarloSimulation.run do
      burn_in 10
      iterations 100
      sample { [1.0, 2.0, 3.0] }
      calculate {|x| [ x[0] * 2.0, x[1] * 2.0, x[2] * 2.0] }
    end    
     
    params = results[0].transpose
    
    nr_params = params.size
    assert_equal nr_params, 3
    
    assert_equal params[0].size, 100
    assert_equal params[1].size, 100
    assert_equal params[2].size, 100
     
    mean = params[0].inject(0.0) {|sum,x| sum + x} / params[0].size
    assert_in_epsilon 2.0, mean, 0.000001   
    
    mean = params[1].inject(0.0) {|sum,x| sum + x} / params[1].size
    assert_in_epsilon 4.0, mean, 0.000001   
    
    mean = params[2].inject(0.0) {|sum,x| sum + x} / params[2].size
    assert_in_epsilon 6.0, mean, 0.000001   
  end
  
  def test_simulation_measurements
    results = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    
    m = Rubybayes::MonteCarloSimulation::Measurement.new(results.shuffle)    
    assert_in_epsilon m.mean, 5.5, 0.000001, "m.mean"
    assert_in_epsilon m.variance, 9.166667, 0.000001, "m.var"
    assert_in_epsilon m.sd, 3.02765, 0.000001, "m.sd"    
    assert_in_epsilon m.quantile(0.75), 7.75 , 0.01 , "m.quantile"      
    
    m2 = Rubybayes::MonteCarloSimulation::Measurement.new(results) 
    
    cum_mean = m2.cumulative_mean
    cum_mean_ref = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
    diff_sum = 0.0
    cum_mean.size.times{|i| diff_sum += (cum_mean[i]-cum_mean_ref[i])**2}
    assert_equal diff_sum, 0.0, "cumulative_mean"
    
    assert_in_epsilon m2.cor(m2), 1.0 , 0.01, "m.acf, lag = 0"
    assert_in_epsilon m2.acf(0), 1.0 , 0.01, "m.acf, lag = 0"
    assert_in_epsilon m2.acf(2), 0.412 , 0.01, "m.acf, lag = 2"
    assert_in_epsilon m2.acf(6), -0.376 , 0.01, "m.acf, lag = 6"
  end
  
end
