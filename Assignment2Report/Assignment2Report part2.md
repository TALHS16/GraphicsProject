## 1. Implement the ability to add point lights to the scene and material to models


![title](images/Ex1.png)

```
now i choise ambiante bordo to the banana, but he is still black,
 its couse the ambient ligth vector is zero's,
 so the final color is zero -> black.
```

![title](images/Ex1_black.png)

```
now i choise ambiante white red to the banana, and pure red to the light,
and we see that the banana is only PURE red,
it is couse the zero's in the light ambient in the green and blue make the final color only red.
```

![title](images/Ex1_red.png)


## 2. Compute diffuse

![title](images/Ex2.png)

## 3. Compute diffuse Gouraud model

![title](images/Gouraud_Ex3.png)

## 4. computes the reflection vector
```
The only light in both images is located at the top right of the frame (about where the yellow circle is).
```

![title](images/Ex4_r.png)

## computes the reflection vector

![title](images/Ex4_l.png)


## 5. Compute the specular light with Phong shading

```
This is only the specular part.
```
![title](images/Ex5_specular.png)

```
And this is the FULL phong model with specular part.
```
![title](images/Ex5_specular2.png)

## 6. Demonstrate different models

![title](images/Ex6_bunny.png)

![title](images/Ex6_cow.png)

![title](images/Ex6_beto.png)