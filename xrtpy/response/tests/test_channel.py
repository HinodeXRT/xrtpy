from xrtpy.response.channel import Channel

def test_creating_a_channel():
    Channel()
    
    
def test_channel_name():
    channel_name = "Al-mesh"
    channel = Channel(channel_name)
    assert channel.name == channel_name