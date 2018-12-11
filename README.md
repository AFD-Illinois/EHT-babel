# EHT-babel
iharm &lt;-> KORAL interface

```bash
./babel path/to/define.h path/to/koral.dat
```

You may need to comment/uncomment the type of KORAL data you expect to read. Do this by modifying the lines in ```babel.c```

```c
geom.dtype = KORAL_GRMHD;
geom.dtype = KORAL_RADGRMHD;
```


