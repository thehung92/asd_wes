# docker file

## build from dockerfile and name the image as 'thehung92phuyen/biotools:v3.0'

```
# cd [path/contain/Dockerfile]
docker build -t thehung92phuyen/biotools:v3.0 .
```

## if change is made to the container and you wish to preserve change
```
docker commit -m 'lib for R visualization' asd_wes thehung92phuyen/biotools:v3.0
docker push thehung92phuyen/biotools:v3.0
```