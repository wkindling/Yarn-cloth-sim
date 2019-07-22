#pragma once
#ifndef CAMERA_H
#define CAMERA_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <iostream>

//Camera Paras
enum Camera_Movement
{
	FORWARD,
	BACKWARD,
	LEFT,
	RIGHT
};

const float YAW = -90.0f;
const float PITCH = 0.0f;
const float SPEED = 2.5f;
const float MOUSE_SENSITIVITY = 0.08f;
const float ZOOM = 48.0f;

//Class Definition
class Camera
{
public:
	Camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 world_up = glm::vec3(0.0f, 0.0f, 1.0f), float yaw = YAW, float pitch = PITCH) :
		Speed(SPEED), MouseSensitivity(MOUSE_SENSITIVITY), Zoom(ZOOM)
	{
		Position = position;
		WorldUp = world_up;
		Yaw = yaw;
		Pitch = pitch;

		setCameraFrame();
	}

	virtual ~Camera() {}

	glm::mat4 getViewMatrix()
	{
		return glm::lookAt(Position, Position + Front, Up);
	}

	void translate(Camera_Movement direction, float delta_t)
	{
		float translation = Speed * delta_t;
		switch (direction)
		{
		case FORWARD: Position += Front * translation;
			break;
		case BACKWARD:Position -= Front * translation;
			break;
		case LEFT:	  Position -= Right * translation;
			break;
		case RIGHT:   Position += Right * translation;
			break;
		default: std::cout << "Camera Translation Information Wrong! " << std::endl;
			break;
		}
	}

	void rotate(float delta_x, float delta_y)
	{
		float delta_yaw = delta_x * MouseSensitivity;
		float delta_pitch = delta_y * MouseSensitivity;

		Yaw += delta_yaw;
		Pitch += delta_pitch;

		if (Pitch > 88.0f) Pitch = 88.0f;
		if (Pitch < -88.0f) Pitch = -88.0f;

		setCameraFrame();
	}

	void zoom(float delta_zoom)
	{
		if (Zoom >= 1.0f && Zoom <= 45.0f) Zoom -= delta_zoom;
		if (Zoom < 1.0f) Zoom = 1.0f;
		if (Zoom > 45.0f) Zoom = 45.0f;
	}


private:
	void setCameraFrame()
	{
		//Pitch y  Yaw z
		glm::vec3 front;
		front.y = cos(glm::radians(Pitch))*cos(glm::radians(Yaw));
		front.z = sin(glm::radians(Pitch));
		front.x = cos(glm::radians(Pitch))*sin(glm::radians(Yaw));

		Front = glm::normalize(front);
		Right = glm::normalize(glm::cross(Front, WorldUp));
		Up = glm::normalize(glm::cross(Right, Front));
	}

public:
	glm::vec3 Position;

	glm::vec3 Front;
	glm::vec3 Up;
	glm::vec3 Right;

	glm::vec3 WorldUp;

	float Yaw;
	float Pitch;

	float Speed;
	float MouseSensitivity;
	float Zoom;

};



#endif