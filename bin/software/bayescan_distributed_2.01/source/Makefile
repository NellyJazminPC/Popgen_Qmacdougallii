# BayeScan makefile

bayescan2: start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 
	g++ -o bayescan2 -static start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 

start.o: start.cpp errors.cpp anyoption.h global_defs.h
	g++ -c start.cpp errors.cpp 

beta.o: beta.cpp global_defs.h
	g++ -c beta.cpp 
      
dirichlet.o: dirichlet.cpp global_defs.h
	g++ -c dirichlet.cpp 

RJupdates.o: RJupdates.cpp global_defs.h
	g++ -c RJupdates.cpp 

MHupdates.o: MHupdates.cpp global_defs.h
	g++ -c MHupdates.cpp 

likelihood.o: likelihood.cpp global_defs.h
	g++ -c likelihood.cpp 

read_write.o: read_write.cpp errors.cpp global_defs.h
	g++ -c read_write.cpp errors.cpp 

anyoption.o: anyoption.cpp anyoption.h 
	g++ -c anyoption.cpp 

clean: 
	rm *.o bayescan2
