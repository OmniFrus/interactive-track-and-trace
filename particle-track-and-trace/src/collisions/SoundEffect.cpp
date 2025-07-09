#include "SoundEffect.h"
#include <iostream>

using namespace std;

SoundEffect::SoundEffect(std::string datapath, const std::shared_ptr<ParticleCollisionCallback> &callback)
    : buffer(), sound(buffer), callback(callback) {
    
    if (!buffer.loadFromFile(datapath)) {
        cerr << "Warning: Sound file not found at \"" << datapath << "\". Sound disabled." << endl;
    }

    sound.setVolume(50);
}

void SoundEffect::handleCollision(int index) {
  callback->handleCollision(index);
  sound.stop();
  sound.play();
}
