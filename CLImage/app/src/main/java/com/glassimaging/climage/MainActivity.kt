package com.glassimaging.climage

import androidx.appcompat.app.AppCompatActivity
import android.os.Bundle
import com.glassimaging.climage.databinding.ActivityMainBinding

class TestCLImage : AppCompatActivity() {

    private lateinit var binding: ActivityMainBinding

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)

        binding = ActivityMainBinding.inflate(layoutInflater)
        setContentView(binding.root)

        // Example of a call to a native method
        binding.sampleText.text = testCLImage()
    }

    /**
     * A native method that is implemented by the 'climage' native library,
     * which is packaged with this application.
     */
    external fun testCLImage(): String

    companion object {
        // Used to load the 'climage' library on application startup.
        init {
            System.loadLibrary("climage")
        }
    }
}