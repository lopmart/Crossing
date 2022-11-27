//#include <GL\freeglut.h>		// Need to include it here because the GL* types are required
//#include "camera.h"
//#define PI 3.1415265359
#define PIdiv180 3.1415265359/180.0

class MousePicker
{
public:
	//private:
	int RECURSION_COUNT = 200;
	float RAY_RANGE = 600;
	glm::vec3 currentRay;
	mat4 projectionMatrix, viewMatrix;
	//CCamera camera;
	glm::vec3 camPos;
	mat4 viewM;

	//public:
	MousePicker() {}

	MousePicker(glm::vec3 pos) {
		//camera = cam;
		camPos = pos;
		projectionMatrix = glm::perspective(60.0f, (GLfloat)glutGet(GLUT_WINDOW_HEIGHT) / (GLfloat)glutGet(GLUT_WINDOW_WIDTH), 1.0f, 256.0f);
		viewMatrix = createViewMatrix(pos);
		viewM = glm::mat4(1.0f);
		viewM = glm::translate(viewM, glm::vec3(0.0f, 0.0f, pos.z));
	}

	mat4 createViewMatrix(glm::vec3 pos)
	{
		viewMatrix = glm::lookAt(
			glm::vec3(pos.x, pos.y, pos.z),
			glm::vec3(0, 0, 0),
			glm::vec3(0, 1, 0));
		viewM = glm::mat4(1.0f);
		viewM = glm::translate(viewM, glm::vec3(0.0f, 0.0f, pos.z));
		return viewMatrix;
	}

	void MousePick(glm::vec3 pos) {
		camPos = pos;
		projectionMatrix = glm::perspective(glm::radians(60.0f), (GLfloat)glutGet(GLUT_WINDOW_HEIGHT) / (GLfloat)glutGet(GLUT_WINDOW_WIDTH), 0.1f, 256.0f);
		viewMatrix = createViewMatrix(pos);
	}

	void update(glm::vec3 pos) {
		//projectionMatrix = glm::perspective(glm::radians(60.0f), (GLfloat)glutGet(GLUT_WINDOW_HEIGHT) / (GLfloat)glutGet(GLUT_WINDOW_WIDTH), 0.1f, 256.0f);
		viewMatrix = createViewMatrix(pos);
		//viewMatrix = viewM;
		currentRay = calculateMouseRay();
	}

	glm::vec3 calculateMouseRay() {
		float mX = mouseX;
		float mY = mouseY;
		glm::vec2 normalizedCoords = getNormalisedDeviceCoordinates(mX, mY);
		glm::vec4 clipCoords = glm::vec4(normalizedCoords.x, normalizedCoords.y, -1.0f, 1.0f);
		glm::vec4 eyeCoords = toEyeCoords(clipCoords.x, clipCoords.y, clipCoords.z, clipCoords.w);
		glm::vec3 worldRay = toWorldCoords(eyeCoords.x, eyeCoords.y, eyeCoords.z, eyeCoords.w);
		return worldRay;
	}

	glm::vec3 toWorldCoords(float x1, float y1, float z1, float w1)
	{
		glm::vec4 eyeCoords = glm::vec4(x1, y1, z1, w1);
		glm::mat4 invertedView = glm::inverse(viewMatrix);
		glm::vec4 transformedVector = invertedView * eyeCoords;
		glm::vec4 rayWorld = glm::vec4(transformedVector.x, transformedVector.y, transformedVector.z, transformedVector.w);
		glm::vec3 mouseRay = glm::vec3(rayWorld.x, rayWorld.y, rayWorld.z);
		mouseRay = glm::normalize(mouseRay);
		return mouseRay;
	}

	glm::vec2 getNormalisedDeviceCoordinates(float mouseX, float mouseY) {
		float x = (2.0f * mouseX) / glutGet(GLUT_WINDOW_WIDTH) - 1.0f;
		float y = (2.0f * mouseY) / glutGet(GLUT_WINDOW_HEIGHT) - 1.0f;
		return glm::vec2(x, y);
	}

	glm::vec4 toEyeCoords(float x1, float y1, float z1, float w1) {
		glm::vec4 clipCoords = glm::vec4(x1, y1, z1, w1);
		glm::mat4 invertedprojection = glm::inverse(projectionMatrix);
		glm::vec4 transformedVector = invertedprojection * clipCoords;
		glm::vec4 eyeCoords = glm::vec4(transformedVector.x, transformedVector.y, transformedVector.z, transformedVector.w);
		return glm::vec4(eyeCoords.x, eyeCoords.y, -1, 0);
	}
};
