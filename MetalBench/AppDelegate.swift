//
//  AppDelegate.swift
//  MetalBench
//
//  Created by Alia on 30/07/2020.
//

import Cocoa
import SwiftUI

// @main
@NSApplicationMain
class AppDelegate: NSObject, NSApplicationDelegate, NSAlertDelegate, NSWindowDelegate {

	var window: NSWindow!
	let renderer = Renderer()

	func applicationDidFinishLaunching(_ aNotification: Notification) {
		
		// Test to ensure a metal device is available...
		if MTLCreateSystemDefaultDevice() == nil {
			// The app has an "oh shit" moment
			let alert = NSAlert.init()
			alert.alertStyle = .critical
			alert.messageText = "No Metal capable GPU available"
			alert.informativeText = "Time to upgrade your Mac?"
			
			alert.addButton(withTitle: "Quit")
			alert.addButton(withTitle: "Fix this error")
			
			let okButton = alert.buttons[0]
			okButton.target = self
			okButton.action = #selector(AppDelegate.killApp)
			let fixButton = alert.buttons[1]
			fixButton.target = self
			fixButton.action = #selector(AppDelegate.fixError)
			
			alert.delegate = self
			
			alert.runModal()
			
		} else {
			// Create the SwiftUI view that provides the window contents.
			let contentView = ContentView()
				.environmentObject(renderer)
			
			// Create the window and set the content view.
			window = NSWindow(
				contentRect: NSRect(x: 0, y: 0, width: 480, height: 300),
				styleMask: [.titled, .closable, .miniaturizable, .fullSizeContentView],
				//			styleMask: [.titled, .closable, .miniaturizable, .resizable, .fullSizeContentView],
				backing: .buffered, defer: false)
			window.titlebarAppearsTransparent = true
			window.isReleasedWhenClosed = false
			window.center()
			window.setFrameAutosaveName("Main Window")
			window.contentView = NSHostingView(rootView: contentView)
			window.makeKeyAndOrderFront(nil)
			window.delegate = self
		}
	}

	func applicationWillTerminate(_ aNotification: Notification) {
		// Insert code here to tear down your application
	}

	
	// Alert button methods
	
	@objc func killApp() {
		exit(EXIT_FAILURE)
	}
	
	@objc func fixError() {
		guard let url = URL(string: "http://www.apple.com/mac/") else {
			exit(EXIT_FAILURE)
		}
		
		NSWorkspace.shared.open(url)
		exit(EXIT_FAILURE)
	}

	// Window delegate methods
	
	func windowWillClose(_ notification: Notification) {
		// Quit on close
		exit(EXIT_SUCCESS)
	}
}

