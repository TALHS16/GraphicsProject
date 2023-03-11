## 1. Orthographic  Projection
### The Original mesh
![title](images/blob.png)

### The Mesh In A Different View Volume
![title](images/blob2.png)

We allowed the user to change the view volume of the camera with the sliders.

## 2. reposition the camera: Camera Frame and World Frame
    The transformation of the camera is exactly like in the meshModel and 
    therefore we created a new class called Transform.
    The calsses MeshModel and Camera inherited from Transform, 
    and Transform has all the transformation (scale rotatioan and translation) 
    of the world frame and the local frame.
    We copied the code from Assignment1B that relate to the meshModel transformation 
    and tranfer it to the new class so that the transformation will be the same 
    for the camera and the MeshModel.
### Video Of The Tranformation


https://user-images.githubusercontent.com/100534541/208247964-f2b5a7c6-ecab-4f36-b788-37e6284f00cf.mp4


## 3. Set The View Volume Automatically 
	we needed the set the view volume automatically when we change the windows size.
    In the main class we checked if the window size has change and if it did we did the following code:
``` cpp
    if (frameBufferWidth != renderer.GetViewportWidth() || frameBufferHeight != renderer.GetViewportHeight())
	{
        // set the width and height in the renderer 
		renderer.SetViewportWidth(frameBufferWidth); 
		renderer.SetViewportHeight(frameBufferHeight);
        //init the renderer buffer color and init the renderer to prevent 
        //from program to collapse when we make the window bigger
		renderer.refreshBuff();
        // set the width and height in the renderer
		Camera::screen_width = frameBufferWidth;
		Camera::screen_hight = frameBufferHeight;
        // set the orth array by the new windows size
		for (int i = 0; i < scene.GetCameraCount(); i++)
		{
			scene.GetCamera(i).refresh();
		}
	} 	
```
### Sanity check: 
    we changed the windows size and checked if the parameters of the view volume changed to.
    and also see if the meshModel centered.
### Original Mesh 
![title](images/Assigment1C_q3_Original.png)
### Width > Height
![title](images/Assigment1C_q3_2.png) 
### Width = Height
![title](images/Assigment1C_q3_3.png)
### Width < Height
![title](images/Assigment1C_q3_4.png) 



## 4. Draw model with world(Green) and local(Red) axes

![title](images/model_with_axes.png)


## 5. Draw model with world(Green) and local(Red) axes
has we can see we rotate the world of the cow and the local axes rotate also as we expact, and the banana axes does't change because we rotate the banana localy.

![title](images/2models_world_local_frame.png)

now we rotate the banana and the cow boht in the world frame so the local axes (in red) also rotated.

![title](images/2models_world_world_frame.png)

## 6. Bounding box

as we can see when we rotate *localy* the object the bouding box change that the bottom line of the box always paralel to the X axie.
![title](images/2bananas_rotate_local.png)


so when we rotate the world of the object the bounding box also rotate and keep paralel to the local X axes of the object.

![title](images/banana_rotated_in_world.png)


## 7. draw normals faces and vertices

for the sanity check we draw a cube and check the normals of the faces:

![title](images/cude_normal_faces.png)

now we draw the normals of the banana object:

![title](images/banana_normals_vertices.png)


## 8. Orthographic Projection VS. The Perspective Projection
    In the video bellow you can see that we implemented perspective projection 
	and the comperison between  the orthographic projection and the perspective projection.
    for the sanity check we check if in the perspective projection we receive a rounded image
	of the mesh that reflect the depth of the mesh.

https://user-images.githubusercontent.com/100534541/208247985-0a2b0ffd-1555-411b-b70d-ed71259b99fb.mp4


## 9. Zoom VS. Getting Closer In Perspective
### Original MeshModel
![title](images/Assigment1C_q9_original.png)
### zoom (by changing the size of the frustum)
![title](images/Assigment1C_q9_Zoom.png)
### Getting Closer
![title](images/Assigment1C_q9_GetingCloser.png)

## 10. LookAt
    We set the camera position to (500,500,500) and orient it towards the origin using 
    LookAt (0,0,1000). 
### The Original Mesh Position
![title](images/Assigment1C_q10_origin.png)

### Video Of The Transformation From (500,500,500) To The Origin (0,0,1000)


https://user-images.githubusercontent.com/100534541/208248024-6a33050a-8137-4490-9945-f54da256816e.mp4


## 11. Dolly Zoom
	As we learned is class the dolly zoom is when the camera appears to simultaneously 
	zoom in and pull back from the subject.
	In the assignment we increased the fovy (the first component in the pers_vec) 
	and decreased the eyeZ to achieve this effect.
https://user-images.githubusercontent.com/100534541/208297014-9d36af35-3fff-4a0d-810f-07097303d9f0.mp4

## 12. Multiple Models & Cameras

https://user-images.githubusercontent.com/100534541/208248045-6da22d8a-b417-4bea-b58f-b06c51fe8634.mp4


## 13. GUI
![title](images/Assigment1C_q12_1.png)
![title](images/Assigment1C_q12_2.png)
![title](images/Assigment1C_q12_3.png)
![title](images/Assigment1C_q12_4.png)
![title](images/Assigment1C_q12_5.png)
![title](images/Assigment1C_q12_6.png)
![title](images/Assigment1C_q12_7.png)
![title](images/Assigment1C_q12_8.png)

Now we see the our full controllers of object and camera:
![title](images/full_controlers.png)

## 14. Additional Features
    we added camera control with the keyboard and added reset button to reset the camera to the initial cordinates.

https://user-images.githubusercontent.com/100534541/208248067-9dd1c8d6-d54f-4a09-973b-271e2fe8be39.mp4
