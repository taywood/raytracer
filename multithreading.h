// * Source from Ryan Westwood (by email) starts here: * //
#pragma once
#include <condition_variable>
#include <functional>
#include <vector>
#include <thread>
#include <queue>

class ThreadPool {
public:
	using Task = std::function<void()>;

	ThreadPool(short threadCount) {
		Start(threadCount);
	}

	~ThreadPool() {
		Stop();
	}

	void Enqueue(Task task) {
		{
			std::unique_lock<std::mutex> lock(mMutex);
			mTasks.emplace(std::move(task));
		}

		mCondition.notify_one();
	}

private:
	std::vector<std::thread> mThreads;
	std::condition_variable mCondition;
	std::mutex mMutex;
	std::queue<Task> mTasks;

	bool mRunning = false;
	void Start(short threadCount) {
		for (short i = 0; i < threadCount; i++)
		{
			mThreads.emplace_back([=] {
				while (true) {
					Task task;
					{
						std::unique_lock<std::mutex> lock(mMutex);

						mCondition.wait(lock, [=] { return mRunning || !mTasks.empty(); });

						if (mRunning && mTasks.empty())
							break;

						task = std::move(mTasks.front());
						mTasks.pop();
					}
					task();
				}
				});
		}
	}

	void Stop() {
		{
			std::unique_lock<std::mutex> lock(mMutex);
			mRunning = true;
		}

		mCondition.notify_all();

		for (auto& thread : mThreads)
			thread.join();
	}
};

// * Source from Ryan Westwood ends here * //

void LineRender(SDL_Surface* screen, hittable_list world, int y, int spp, int max_depth, camera* cam) {
	const float aspect_ratio = 16.0 / 9;
	const int image_width = screen->w;
	const int image_height = static_cast<int>(image_width / aspect_ratio);

	const Colour black(0, 0, 0);
	Colour pix_col(black);

	for (int x = 0; x < screen->w; ++x) {
		pix_col = black;
		for (int s = 0; s < spp; s++) {
			auto u = double(x + random_double()) / (image_width - 1);
			auto v = double(y + random_double()) / (image_height - 1);
			ray ray = cam->get_ray(u, v);
			pix_col = pix_col + ray_colour(ray, world, max_depth);
		}
		pix_col /= 255.f * spp;
		pix_col.x = sqrt(pix_col.x);
		pix_col.y = sqrt(pix_col.y);
		pix_col.z = sqrt(pix_col.z);
		pix_col *= 255;
		Uint32 colour = SDL_MapRGB(screen->format, pix_col.x, pix_col.y, pix_col.z);
		putpixel(screen, x, y, colour);
	}
}