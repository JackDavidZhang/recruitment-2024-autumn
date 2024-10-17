### init

#### Local WSL

- Tiny 258.552 ms
- Small 1.032798 s

#### Remote RaspberryPi4B

- Tiny 19.070420 s

- Small 81.802380 s

### O3

#### Remote RaspberryPi4B

- Tiny 1.360060 s
- Small 5.503337 s

### Compress DP Matrix

#### Local WSL

- Tiny 239.759000 ms
- Small 963.333000 ms
- Medium 98.535538 s

#### Remote RaspberryPi4B

- Tiny 1.258235 s
- Small 5.079680 s
- Medium 633.32 s

### Multithread(std::thread)

#### Local WSL

- Tiny 41.169000 ms
- Small 163.61600 ms
- Medium 16.926868 s

#### Remote RaspberryPi4B

- Tiny 361.768000 ms
- Small 1.489942 s
- Medium 705.848305 s

##### Single Thread(Medium)

- 279.216230 s
- 282.740769 s
- 282.880155 s
- 285.955705 s
- 275.629637 s
- 281.382652 s
- 287.613684 s
- 288.292549 s
- 112.112895 s
- 147.349991 s

### Multithread(OpenMP)

#### Local WSL

- Tiny 48.554000 ms
- Small 164.031000 ms
- Medium 16.965020 s

#### Remote RaspberryPi4B

- Tiny 357.441000 ms
- Small 1.486954 s
- Medium 779.295774 s

### Multithread(OpenMP SIMD)

#### Local WSL

- Tiny 50.437000 ms
- Small 143.368000 ms
- Medium 12.877843 s

#### Remote RaspberryPi4B

- Tiny 387.141000 ms
- Small 1.603563 s
- Medium 752.945457 s

#### Remote Server ?

- Tiny 66.345000 ms
- Small 241.631000 ms
- Medium 24.396443 s

### Decouple computing(OpenMP SIMD)

#### Local WSL

- Tiny 59.594000 ms
- Small 185.129000 ms
- Medium 17.385655 s

### Decouple computing(OpenMP)

#### Local WSL

- Tiny 43.637000 ms
- Small 128.511000 ms
- Medium 12.364470 s