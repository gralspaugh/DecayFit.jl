using DecayFit 

irf = simulate_irf(0.047,256,44.7,100000,0.25,1.5)

dk = simulate_dk(0.097,256,44.7,irf,zeros(Int,size(irf)),
[0.1, 2.3, 0, 0, 0, 0],1,0)
