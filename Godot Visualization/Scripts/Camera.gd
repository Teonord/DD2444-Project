extends Node3D

@onready var camera = $Camera3D
@export var sensitivity = 0.0025
@export var speed = 3
	
func _input(event):
	if event.is_action_pressed("m1"):
		Input.mouse_mode = Input.MOUSE_MODE_CAPTURED
	elif event.is_action_pressed("ui_cancel"):
		Input.mouse_mode = Input.MOUSE_MODE_VISIBLE
	elif event is InputEventMouseMotion and Input.MOUSE_MODE_CAPTURED:
		rotate_y(-event.relative.x * sensitivity)
		camera.rotate_x(-event.relative.y * 0.01)
		camera.rotation.x = clamp(camera.rotation.x, -PI/2, PI/2)
	
func _process(delta):
	var indir = Input.get_vector("left", "right", "forward", "backward")
	var vert = Input.get_axis("down", "up")
	position += transform.basis * Vector3(indir.x, 0, indir.y).normalized() * speed * delta
	position.y += vert * speed * delta
	
