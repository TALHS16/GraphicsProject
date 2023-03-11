## 1. Draw the bounding rectangle
```
For drawing the bounding rectangle we need to enable that option.
We calculate the color of the ractangle in relation to the depth.
The color black represents it's closeness to the screen and
the color white represents it's farness to the screen.
In every triangle we calculate the avarage of the z cordinate and
check how far it is from the min z cordinate.
```

![title](images/img3.png)

## 2. Pass
```
we chose to implement the edge walking approach.
```

## 3. Fill the traingle with random color
```
To fill the traingles with color we enable the 'Fill Traingle' option.
```
![title](images/img4.png)


## 4. Z - buffer
```
First we can see the z buffer algoritem on one model:
```
![title](images/img1.png)
```
It is possible to recognize that the back of the cow is darker than the front as we expected to see.
```
### Bunny model
```
The bunny model using color buffer:
```
![title](images/img4.png)
```
The bunny model using z buffer:
```
![title](images/img5.png)

### Camera model

```
The camera model using color buffer:
```
![camera_color_buffer](https://user-images.githubusercontent.com/100534541/211156255-97d6719b-1a70-43b1-bc4c-0f7220eb355a.png)

```
The camera model using z buffer:
```
![camera_z_buffer](https://user-images.githubusercontent.com/100534541/211156262-04d94f4c-df7c-4b5a-8f56-7534d96328f8.png)

### Banana model

```
The Banana model using color buffer:
```
![banana_color_buffer](https://user-images.githubusercontent.com/100534541/211156281-ecf28d05-b0f8-408d-96e7-aa224ee6bd79.png)

```
The Banana model using z buffer:
```
![banana_z_buffer](https://user-images.githubusercontent.com/100534541/211156291-6dcd5fae-eadf-454f-8d02-082edfb5d503.png)
