#define _USE_MATH_DEFINES
#include <cmath>
#include <imgui/imgui.h>
#include <stdio.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <nfd.h>
#include <iostream>

#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Renderer.h"
#include "Scene.h"
#include "Utils.h"
#include "MeshModel.h"

/**
 * Fields
 */
bool show_z_buffer = false;
bool show_another_window = false;
bool flag_model = true;
glm::vec4 clear_color = glm::vec4(0.8f, 0.8f, 0.8f, 1.00f);
glm::vec3 model_color = glm::vec3(0.0f, 0.0f, 0.0f);
int windowWidth = 1280, windowHeight = 720;


/**
 * Function declarations
 */
static void GlfwErrorCallback(int error, const char* description);
GLFWwindow* SetupGlfwWindow(int w, int h, const char* window_name);
ImGuiIO& SetupDearImgui(GLFWwindow* window);
void StartFrame();
void RenderFrame(GLFWwindow* window, Scene& scene, Renderer& renderer, ImGuiIO& io);
void Cleanup(GLFWwindow* window);
void DrawImguiMenus(ImGuiIO& io, Scene& scene);
int frameBufferWidth, frameBufferHeight;
/**
 * Function implementation
 */
void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset)
{
	ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
	// TODO: Handle mouse scroll here
}

int main(int argc, char **argv)
{
    
	GLFWwindow* window = SetupGlfwWindow(windowWidth, windowHeight, "Mesh Viewer");
    GLFWwindow* z_buffer_window=NULL;
	if (!window)
		return 1;
	
	glfwMakeContextCurrent(window);
	glfwGetFramebufferSize(window, &frameBufferWidth, &frameBufferHeight);
	Renderer renderer = Renderer(frameBufferWidth, frameBufferHeight);
	Scene scene = Scene();
	std::shared_ptr<Camera> firstCamera (new Camera());
	Camera::screen_width = frameBufferWidth;
	Camera::screen_hight = frameBufferHeight;
	firstCamera->SetCameraLookAt();
	firstCamera->SetProjection();
	scene.AddCamera(firstCamera);
	scene.SetActiveCameraIndex(0);

	//  std::shared_ptr<MeshModel> firstModel = Utils::LoadMeshModel("/Users/idoshitrit/Desktop/computer-graphics-2023-ido-tal/Data/cow.obj");
	//  scene.AddModel(firstModel);
	//   std::shared_ptr<MeshModel> secondtModel = Utils::LoadMeshModel("/Users/idoshitrit/Desktop/computer-graphics-2023-ido-tal/Data/banana.obj");
	//   scene.AddModel(secondtModel);
	
	ImGuiIO& io = SetupDearImgui(window);
	glfwSetScrollCallback(window, ScrollCallback);
    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
		StartFrame();
		DrawImguiMenus(io, scene);
		RenderFrame(window, scene, renderer, io);
        renderer.get_show_z_buffer() = show_z_buffer;
//        if (show_z_buffer){
////            if(z_buffer_window==NULL){
////                z_buffer_window = SetupGlfwWindow(windowWidth, windowHeight, "Z Buffer");
////                glfwMakeContextCurrent(z_buffer_window);
////                glfwGetFramebufferSize(z_buffer_window, &frameBufferWidth, &frameBufferHeight);
////
////            }
////            renderer.draw_Z_buffer(z_buffer_window);
////            glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, windowWidth, windowHeight,
////                         0, GL_RED, GL_UNSIGNED_BYTE, renderer.zBuffer);
//
//        }
//        else{
////            if(z_buffer_window!=NULL){
////                glfwDestroyWindow(z_buffer_window);
////            }
//            z_buffer_window=NULL;
//
//        }
    }

	Cleanup(window);
    return 0;
}

static void GlfwErrorCallback(int error, const char* description)
{
	fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

GLFWwindow* SetupGlfwWindow(int w, int h, const char* window_name)
{
	glfwSetErrorCallback(GlfwErrorCallback);
	if (!glfwInit())
		return NULL;
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	
	#if __APPLE__
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	#endif
	
	GLFWwindow* window = glfwCreateWindow(w, h, window_name, NULL, NULL);
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync
						 // very importent!! initialization of glad
						 // https://stackoverflow.com/questions/48582444/imgui-with-the-glad-opengl-loader-throws-segmentation-fault-core-dumped

	gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
	return window;
}

ImGuiIO& SetupDearImgui(GLFWwindow* window)
{
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init();
	ImGui::StyleColorsDark();
	return io;
}

void StartFrame()
{
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();
}

void RenderFrame(GLFWwindow* window, Scene& scene, Renderer& renderer, ImGuiIO& io)
{
	ImGui::Render();
	int frameBufferWidth, frameBufferHeight;
	glfwMakeContextCurrent(window);
	glfwGetFramebufferSize(window, &frameBufferWidth, &frameBufferHeight);
	if (frameBufferWidth != renderer.GetViewportWidth() || frameBufferHeight != renderer.GetViewportHeight())
	{

		renderer.SetViewportWidth(frameBufferWidth);
		renderer.SetViewportHeight(frameBufferHeight);
		renderer.refreshBuff();
		Camera::screen_width = frameBufferWidth;
		Camera::screen_hight = frameBufferHeight;
		for (int i = 0; i < scene.GetCameraCount(); i++)
		{
			scene.GetCamera(i).refresh();
		}

		// TODO: Set new aspect ratio
	}

	if (!io.WantCaptureKeyboard)
	{
		if(scene.GetModelCount()>0){
			MeshModel& activeModel = scene.GetActiveModel();
			// TODO: Handle keyboard events here
			if (io.KeysDown[io.KeyMap[ImGuiKey_UpArrow]]) // 'U'
			{
				activeModel.get_move_factor()[1]++;
			}
			if (io.KeysDown[io.KeyMap[ImGuiKey_DownArrow]]) // 'D'
			{
				activeModel.get_move_factor()[1]--;
			}
			if (io.KeysDown[io.KeyMap[ImGuiKey_LeftArrow]]) // 'L'
			{
				activeModel.get_move_factor()[0]--;
			}
			if (io.KeysDown[io.KeyMap[ImGuiKey_RightArrow]]) // 'R'
			{
				activeModel.get_move_factor()[0]++;
			}

			if (io.KeysDown[87]|| io.KeysDown[119]) // 'W'
			{
				scene.GetActiveCamera().get_move_factor()[1]+=10;
			}
			if (io.KeysDown[83] || io.KeysDown[115]) // 'S'
			{
				scene.GetActiveCamera().get_move_factor()[1]-=10;
			}
			if (io.KeysDown[68] || io.KeysDown[115]) // 'D'
			{
				scene.GetActiveCamera().get_move_factor()[0]+=10;
			}
			if (io.KeysDown[65] || io.KeysDown[97]) // 'A'
			{
				scene.GetActiveCamera().get_move_factor()[0]-=10;
			}
		}

	}

	if (!io.WantCaptureMouse)
	{
        if(scene.GetModelCount()>0){
            // TODO: Handle mouse events here
            if (io.MouseDown[0])
            {
                
                scene.GetActiveModel().get_scale_factor()++;
                
                // Left mouse button is down
            }
            else if (io.MouseDown[1])
            {
                if(scene.GetActiveModel().get_scale_factor() >= 0)
                    scene.GetActiveModel().get_scale_factor()--;
                // Right mouse button is down
            }
        }
	}

	renderer.ClearColorBuffer(clear_color);
	renderer.Render(scene);
	renderer.SwapBuffers();

	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
	glfwMakeContextCurrent(window);
	glfwSwapBuffers(window);
}

void Cleanup(GLFWwindow* window)
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();
}

void DrawImguiMenus(ImGuiIO& io, Scene& scene)
{
	/**
	 * MeshViewer menu
	 */
	
	ImGui::Begin("MeshViewer Menu");
	
	// Menu Bar
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Open", "CTRL+O"))
			{
				nfdchar_t* outPath = NULL;
				nfdresult_t result = NFD_OpenDialog("obj;", NULL, &outPath);
				if (result == NFD_OKAY)
				{
					//MeshModel model = 
					scene.AddModel(Utils::LoadMeshModel(outPath));

					free(outPath);
				}
				else if (result == NFD_CANCEL)
				{
				}
				else
				{
				}

			}
			ImGui::EndMenu();
		}
		// show window for each MeshModel
		// loop over all the models

		bool focus;

		for(int i=0; i<scene.GetModelCount();i++){
			MeshModel& model = scene.GetModel(i);
			ImGui::Begin(("Control panel to MeshModel " + std::to_string(i) + ":" + model.GetModelName()).c_str(), &flag_model);
			
			ImGui::ColorEdit3("Change color", (float*)&model.get_model_color());
			focus = ImGui::IsItemActivated();
            
            ImGui::Checkbox("Show Object", &model.get_show_object());
            focus = focus || ImGui::IsItemActivated();
            
            ImGui::Checkbox("Fill Triangle", &model.get_fill_triangle());
            focus = focus || ImGui::IsItemActivated();

			ImGui::SliderFloat("Scale", &model.get_scale_factor(), 0.0f, 5000.0f);
			focus = focus || ImGui::IsItemActivated();
            
            if (ImGui::CollapsingHeader("Local"))
            {
                ImGui::SliderFloat("Move X", &model.get_move_factor()[0], -windowWidth / 2, windowWidth/2);
                focus = focus || ImGui::IsItemActivated();
                ImGui::SliderFloat("Move Y", &model.get_move_factor()[1], -windowHeight / 2, windowHeight/2);
                focus = focus || ImGui::IsItemActivated();
                //ImGui::SliderFloat("Move Z", &model.get_move_factor()[2], -windowHeight / 2, windowHeight/2);
                ImGui::SliderFloat("Move Z", &model.get_move_factor()[2], -5000.0f, 5000.0f);
                focus = focus || ImGui::IsItemActivated();

                ImGui::SliderFloat("Self Rotate X", &model.get_rotate_factor()[0], 0.0f, 360.0f);
                focus = focus || ImGui::IsItemActivated();

                ImGui::SliderFloat("Self Rotate Y", &model.get_rotate_factor()[1], 0.0f, 360.0f);
                focus = focus || ImGui::IsItemActivated();

                
                ImGui::SliderFloat("Self Rotate Z", &model.get_rotate_factor()[2], 0.0f, 360.0f);
                focus = focus || ImGui::IsItemActivated();
            }

            if (ImGui::CollapsingHeader("World"))
            {
                ImGui::SliderFloat("World Rotate X", &model.get_rotate_world_factor()[0], 0.0f, 360.0f);
                focus = focus || ImGui::IsItemActivated();

                ImGui::SliderFloat("World Rotate Y", &model.get_rotate_world_factor()[1], 0.0f, 360.0f);
                focus = focus || ImGui::IsItemActivated();

                ImGui::SliderFloat("World Rotate Z", &model.get_rotate_world_factor()[2], 0.0f, 360.0f);
                focus = focus || ImGui::IsItemActivated();
                
                ImGui::SliderFloat("World Move X", &model.get_world_move_factor()[0], -windowWidth / 2, windowWidth/2);
                focus = focus || ImGui::IsItemActivated();
                
                ImGui::SliderFloat("World Move Y", &model.get_world_move_factor()[1], -windowWidth / 2, windowWidth/2);
                focus = focus || ImGui::IsItemActivated();
                
                ImGui::SliderFloat("World Move Z", &model.get_world_move_factor()[2], -windowWidth / 2, windowWidth/2);
                focus = focus || ImGui::IsItemActivated();
            }
            
            ImGui::Checkbox("Show Axes", &model.get_show_axes());
            focus = focus || ImGui::IsItemActivated();
            
            ImGui::Checkbox("Show Faces Normals", &model.get_show_faces_normals());
            focus = focus || ImGui::IsItemActivated();
            
            ImGui::Checkbox("Show Vertices Normals", &model.get_show_verteces_normals());
            focus = focus || ImGui::IsItemActivated();
            
            ImGui::Checkbox("Show Triangle Bound", &model.get_show_bound_triangle());
            focus = focus || ImGui::IsItemActivated();

            
            if (ImGui::Button("Reset Model")) // Buttons return true when clicked (most widgets return true when edited/activated)
            {
                model.init_trans();
            }
            
            if (ImGui::Button("Remove Model")) // Buttons return true when clicked (most widgets return true when edited/activated)
            {
                scene.RemoveModelByIndex(i);
            }
			
            
			if(focus){
				//std::cout<<model.GetModelName()<<std::endl;
				scene.SetActiveModelIndex(i);
			}

			ImGui::End();
		}
		// TODO: Add more menubar items (if you want to)
		ImGui::EndMainMenuBar();
	}

	// Controls
	ImGui::ColorEdit3("Clear Color", (float*)&clear_color);
	// TODO: Add more controls as needed
	
	ImGui::End();

	/**
	 * Imgui demo - you can remove it once you are familiar with imgui
	 */
	
	// 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).

	bool flag_state = 0;
	static int camera_counter = 1;
	ImGui::Begin("Camera Controller");                          // Create a window called "Hello, world!" and append into it.

	ImGui::Text("camera state:");

	//ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
	//ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color
	if (ImGui::Button("add camera")) // Buttons return true when clicked (most widgets return true when edited/activated)
	{
		std::shared_ptr<Camera> firstCamera(new Camera());
		firstCamera->SetCameraLookAt();
		firstCamera->SetProjection();
		scene.AddCamera(firstCamera);
		scene.SetActiveCameraIndex(camera_counter);
		camera_counter++;
	}
	ImGui::Text("counter = %d", camera_counter);
	ImGui::Text("current camera index = %d", scene.GetActiveCameraIndex());

	static const char* items_pro[] = { "orthographic", "perspective" };
	int current_item_pro = scene.GetActiveCamera().get_camera_view();

	if (ImGui::Button("prev camera")) // Buttons return true when clicked (most widgets return true when edited/activated)
	{
		int curr = scene.GetActiveCameraIndex();
    scene.SetActiveCameraIndex((curr+scene.GetCameraCount()-1) % scene.GetCameraCount());
	}

	if (ImGui::Button("next camera")) // Buttons return true when clicked (most widgets return true when edited/activated)
	{
		int curr = scene.GetActiveCameraIndex();
        scene.SetActiveCameraIndex((curr+1) % scene.GetCameraCount());
		current_item_pro = scene.GetActiveCamera().get_camera_view();
	}

	ImGui::Combo("Select camera view state", &current_item_pro, items_pro, IM_ARRAYSIZE(items_pro));
	scene.GetActiveCamera().set_camera_view(current_item_pro);

	if (ImGui::CollapsingHeader("Projection"))
	{
		if (current_item_pro == 0)
		{
			ImGui::SliderFloat("Left", &scene.GetActiveCamera().get_ortho()[0], -windowWidth / 2, 0.0f);
			ImGui::SliderFloat("Right", &scene.GetActiveCamera().get_ortho()[1], 0.0f, windowWidth / 2);
			ImGui::SliderFloat("Bottom", &scene.GetActiveCamera().get_ortho()[2], -windowHeight / 2, 0.0f);
			ImGui::SliderFloat("Top", &scene.GetActiveCamera().get_ortho()[3], 0.0f, windowHeight / 2);
			ImGui::SliderFloat("Near", &scene.GetActiveCamera().get_ortho()[4], 0.0f, windowWidth);
			ImGui::SliderFloat("Far", &scene.GetActiveCamera().get_ortho()[5], 0.0f, windowWidth);
		}
		if (current_item_pro == 1)
		{
			float dolly = scene.GetActiveCamera().get_dolly_zoom();
			ImGui::SliderFloat("Fovy", &scene.GetActiveCamera().get_pers()[0], 0.f, 360.0f);
			ImGui::InputFloat("Aspect", &scene.GetActiveCamera().get_pers()[1], 0.1f, 1.0f, "%.1f");
			ImGui::InputFloat("Near", &scene.GetActiveCamera().get_pers()[2], 1.0f, 1.0f, "%.1f");
			ImGui::InputFloat("Far", &scene.GetActiveCamera().get_pers()[3], 1.0f, 1.0f, "%.1f");
			ImGui::InputFloat("Dolly Zoom", &dolly, 0.1f, 0.1f, "%.1f");
			if (dolly != scene.GetActiveCamera().get_dolly_zoom())
				scene.GetActiveCamera().setDollyZoom(dolly);
		}
	}

	if (ImGui::CollapsingHeader("Transformation"))
	{
		ImGui::SliderFloat("Scale", &scene.GetActiveCamera().get_scale_factor(), 0.0f, 50.0f);

		ImGui::SliderFloat("Move X", &scene.GetActiveCamera().get_move_factor()[0], -windowWidth/2, windowWidth/2);

		ImGui::SliderFloat("Move Y", &scene.GetActiveCamera().get_move_factor()[1], -windowHeight/2, windowHeight/2);

		ImGui::SliderFloat("Move Z", &scene.GetActiveCamera().get_move_factor()[2], -windowHeight/2 , windowHeight/2);

		ImGui::SliderFloat("Self Rotate X", &scene.GetActiveCamera().get_rotate_factor()[0], 0.0f, 360.0f);

		ImGui::SliderFloat("Self Rotate Y", &scene.GetActiveCamera().get_rotate_factor()[1], 0.0f, 360.0f);

		ImGui::SliderFloat("Self Rotate Z", &scene.GetActiveCamera().get_rotate_factor()[2], 0.0f, 360.0f);

		ImGui::SliderFloat("World Rotate X", &scene.GetActiveCamera().get_rotate_world_factor()[0], 0.0f, 360.0f);

		ImGui::SliderFloat("World Rotate Y", &scene.GetActiveCamera().get_rotate_world_factor()[1], 0.0f, 360.0f);

		ImGui::SliderFloat("World Rotate Z", &scene.GetActiveCamera().get_rotate_world_factor()[2], 0.0f, 360.0f);

		ImGui::SliderFloat("World Move X", &scene.GetActiveCamera().get_world_move_factor()[0], -windowWidth/2, windowWidth/2);

		ImGui::SliderFloat("World Move Y", &scene.GetActiveCamera().get_world_move_factor()[1], -windowHeight / 2, windowHeight/2);

		ImGui::SliderFloat("World Move Z", &scene.GetActiveCamera().get_world_move_factor()[2], -windowHeight / 2, windowHeight / 2);
	}
	if (ImGui::CollapsingHeader("Look At"))
	{
		ImGui::InputFloat("Eye X", &scene.GetActiveCamera().get_eye()[0], 1.0f, 1.0f, "%.1f");
		ImGui::InputFloat("Eye Y", &scene.GetActiveCamera().get_eye()[1], 1.0f, 1.0f, "%.1f");
		ImGui::InputFloat("Eye Z", &scene.GetActiveCamera().get_eye()[2], 1.0f, 1.0f, "%.1f");

		ImGui::InputFloat("At X", &scene.GetActiveCamera().get_at()[0], 1.0f, 1.0f, "%.1f");
		ImGui::InputFloat("At Y", &scene.GetActiveCamera().get_at()[1], 1.0f, 1.0f, "%.1f");
		ImGui::InputFloat("At Z", &scene.GetActiveCamera().get_at()[2], 1.0f, 1.0f, "%.1f");

		ImGui::InputFloat("Up X", &scene.GetActiveCamera().get_up()[0], 1.0f, 1.0f, "%.1f");
		ImGui::InputFloat("Up Y", &scene.GetActiveCamera().get_up()[1], 1.0f, 1.0f, "%.1f");
		ImGui::InputFloat("Up Z", &scene.GetActiveCamera().get_up()[2], 1.0f, 1.0f, "%.1f");
	}
	if (ImGui::Button("Reset"))
	{
		scene.GetActiveCamera().init();
	}
	//ImGui::InputFloat("Dolly Zoom", &scene.GetActiveCamera().get_dolly_zoom(), 0.1f, 0.1f, "%.1f");
	
	ImGui::End();

	// 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
	{
		static float f = 0.0f;
		static int counter = 0;

		ImGui::Begin("MeshViewer Menu");                          // Create a window called "Hello, world!" and append into it.

		ImGui::Checkbox("Show Z Buffer", &show_z_buffer);      // Edit bools storing our window open/close state
		//ImGui::Checkbox("Another Window", &show_another_window);

		//ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
		//ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color

//		if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
//			counter++;
//		ImGui::SameLine();
//		ImGui::Text("counter = %d", counter);

		ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
		ImGui::End();
	}

	// 3. Show another simple window.
	if (show_another_window)
	{
		ImGui::Begin("Another Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
		ImGui::Text("Hello from another window!");
		if (ImGui::Button("Close Me"))
			show_another_window = false;
		ImGui::End();
	}

}
