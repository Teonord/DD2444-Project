extends Node3D

var particle_scene = preload("res://Scenes/Bird.tscn")
@onready var json = JSON.new()

func _ready():
	var file = FileAccess.open("res://res", FileAccess.READ)
	var settings = []
	for item in file.get_line().split(" "):
		settings.append(float(item))
	print(settings)
		
	var line = file.get_line()
	var particles = []
	var particle_line = json_parse(line)
	for i in range(settings[0]):
		var particle = particle_scene.instantiate()
		add_child(particle)
		var pos = Vector3(particle_line[i][0], particle_line[i][1], particle_line[i][2])
		particle.position = pos
		particle.look_at(pos + Vector3(particle_line[i][3], particle_line[i][4], particle_line[i][5]))
		particles.append(particle)
		
	for i in range(settings[1] - 1):
		line = file.get_line()
		particle_line = json_parse(line)
		for j in range(settings[0]):
			var pos = Vector3(particle_line[j][0], particle_line[j][1], particle_line[j][2])
			particles[j].position = pos
			particles[j].look_at(pos + Vector3(particle_line[j][3], particle_line[j][4], particle_line[j][5]))
		await get_tree().create_timer(settings[3]).timeout
		
		
func json_parse(line):
	var result = json.parse("[" + line + "]")
	if result == OK:
		return json.data
	else:
		print("Error in parsing!")
		return []
