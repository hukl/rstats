#          Attack     Defense
# Arsenal  0.00000000 -0.29966558
# Wolves  -0.94049578 -0.56096282
# 
# 
# Home Team Winning Chance
# exp(Home_Attack - (Away_Defense) + HomeWinBonus)
# 
# Away Team Winning Chance
# exp(Away_Attack  - (Home_Defense))


SEGMENTS = 1.0

def probabilities
  # home_team  = [ 0.0,        -0.29966558]
  # away_team  = [-0.94049578, -0.56096282]
  home_team  = [ 1.0, 1.0]
  away_team  = [ -1.0, -0.50600000]
  home_bonus = 0.4579831
  
  goal_prob_home = Math.exp(home_team[0] - (away_team[1]) + home_bonus)
  goal_prob_away = Math.exp(away_team[0] - (home_team[1]))
  
  [goal_prob_home/SEGMENTS, goal_prob_away/SEGMENTS]
end


def rpois(probability)
	mu = probability
	m  = [1, mu.to_i].max
	q  = p0 = p = Math.exp(-mu)
	pp = []
	l  = 0

	while true
	  # roll the dice
		u  = rand
		
		# if the random number is lower than p0 than return 0 as a short circuit
		if u <= p0 
			return 0
		end
		
		# if we went through all 35 iterations without finding an answer
		# we're 
		if l != 0
			(u <= 0.458 ? 1 : [l,m].min).upto(l) do |t|
				if pp[t] && u <= pp[t]
					return t
				end
			end	
		end
	
		l+1.upto(35) do |k|
		  p *= mu / k
		  q += p
			pp[k] = q
			if u <= q
			  return k
		  end
		end
		
		l = 35
	end
end


probs   = probabilities
results = []

1000000.times do
  result = []

  SEGMENTS.to_i.times do |i|
    #result[i] = rpois(probs[0])
    result[i] = rpois(20.0)
  end

  results << result.inject(0) {|acc, e| acc += e}
end

File.open("/tmp/sample", "w") do |file|
  file.puts results
end

# 10.times do
#   probs = probabilities
# 
#   result = []
# 
#   SEGMENTS.to_i.times do |i|
#     result[i] = [rpois(probs[0]), rpois(probs[1])]
#   end
# 
# 
#   result.inject([0,0]) { |acc, e|
#     acc[0] += e[0]
#     acc[1] += e[1]
#     acc
#   }
#   
#   puts result.inspect
# end