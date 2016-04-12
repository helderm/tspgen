# tspgen
Genetic parallel algorithm for approximating the Travelling Salesman Problem

### Contributing

* Install lib `mpich2`

```bash
$ sudo apt-get install mpich2
```

* Compile

```bash
$ mpicc -o tspgen *.c -DDEBUG
```

* Run

```bash
$ mpiexec -n 2 ./tspgen <population size> <mutation rate> <num generations> <num elitism> <mutation size> <num breeders>
```