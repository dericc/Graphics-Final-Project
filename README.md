<HTML>
<HEAD>
<BODY>
<CENTER>
<H2>
COS 426, Spring 2014
</H2> 
<H3>
Deric Cheng<BR>       
Matt Goldsmith<BR> 
Riley Thomasson<BR>   
Andrew Werner<BR>  
</H3>                     
</CENTER>
<HR><BR>


<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 
<H2>List of Implemented Features</H2>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 

For our final project, we created a 2.5 dimension platforming game that has many 
of the necessary features for platforming game: coins, moving platforms, enemies that
bounce and track the player, and decorative fire textures (the most important part of making
any video game look sick). The game is built on a .scn file that takes in the parameters for nearly everything in a level, from the skybox to the textures to the coin location, and the game transitions smoothly between levels (e.g. loads a new scene file immediately once the goal is reached). Our most useful and creative feature, though, is a level editor that
allows you to generate platforms by mouse, move platforms, create coins, create enemies, and so forth, and then immediately play the game in the window. <BR>

Here are the features that we implemented, for a total of 4 basic features and 8 extra features
<ul>
<li><b> 3D perspective viewing and objects.</b>
<li><b> Lighting and smooth shading.</b>
<li><b> User Input</b>
  <li><b> Computer Control over some element in the scene</b>
<li><b> 1) Texture Mapping</b>
<li><b> 2) Multiple Views</b>
<li><b> 3) On-Screen Control Panel</b>
<li><b> 4) Procedural and physically-based modeling.</b>  

<li><b> 5) Collision Detection</b>
<li><b> 6) Simulated Dynamics</b>
<li><b> 7) Sound </b>
<li><b> 8) Game Level Editor</b>

</ul>

There are also many features of the game that do not fit under these 8 categories. We list them following our overview of the important features. 

<TABLE>
<TBODY>
  <TR>
    <TD vAlign=top align=middle>
       <A href="screenshots/image2.jpg"><IMG  width=256 src="screenshots/image2.jpg"></A><BR>       
       <A href="screenshots/image2.jpg">screenshots/image2.jpg</A>
    </TD>
    <TD vAlign=top align=middle>
       <A href="screenshots/image3.jpg"><IMG width=256  src="screenshots/image3.jpg"></A><BR>       
       <A href="screenshots/image3.jpg">screenshots/image3.jpg</A>
    </TD>  
    <TD vAlign=top align=middle>
       <A href="screenshots/image4.jpg"><IMG width=256  src="screenshots/image4.jpg"></A><BR>       
       <A href="screenshots/image4.jpg">screenshots/image4.jpg</A>
    </TD>  
        <TD vAlign=top align=middle>
       <A href="screenshots/image5.jpg"><IMG width=256  src="screenshots/image5.jpg"></A><BR>       
       <A href="screenshots/image5.jpg">screenshots/image5.jpg</A>
  </TD>
  </TR>
  <TR>
        <TD vAlign=top align=middle>
       <A href="screenshots/image6.jpg"><IMG width=256  src="screenshots/image6.jpg"></A><BR>       
       <A href="screenshots/image6.jpg">screenshots/image6.jpg</A>
    </TD>  
        <TD vAlign=top align=middle>
       <A href="screenshots/image7.jpg"><IMG width=256  src="screenshots/image7.jpg"></A><BR>       
       <A href="screenshots/image7.jpg">screenshots/image7.jpg</A>
    </TD>  
        <TD vAlign=top align=middle>
       <A href="screenshots/image8.jpg"><IMG width=256  src="screenshots/image8.jpg"></A><BR>       
       <A href="screenshots/image8.jpg">screenshots/image8.jpg</A>
    </TD>  
    <TD vAlign=top align=middle>
       <A href="screenshots/image9.jpg"><IMG width=256  src="screenshots/image9.jpg"></A><BR>       
       <A href="screenshots/image9.jpg">screenshots/image8.jpg</A>
    </TD>  
</TR>
  <TR>
        <TD vAlign=top align=middle>
       <A href="screenshots/image10.jpg"><IMG width=256  src="screenshots/image10.jpg"></A><BR>       
       <A href="screenshots/image10.jpg">screenshots/image6.jpg</A>
    </TD>  
        <TD vAlign=top align=middle>
       <A href="screenshots/image11.jpg"><IMG width=256  src="screenshots/image11.jpg"></A><BR>       
       <A href="screenshots/image11.jpg">screenshots/image7.jpg</A>
    </TD>  
        <TD vAlign=top align=middle>
       <A href="screenshots/image12.jpg"><IMG width=256  src="screenshots/image12.jpg"></A><BR>       
       <A href="screenshots/image12.jpg">screenshots/image8.jpg</A>
    </TD>  
    <TD vAlign=top align=middle>
       <A href="screenshots/image13.jpg"><IMG width=256  src="screenshots/image13.jpg"></A><BR>       
       <A href="screenshots/image13.jpg">screenshots/image8.jpg</A>
    </TD>  
</TR>
  <TR>
        <TD vAlign=top align=middle>
       <A href="screenshots/image14.jpg"><IMG width=256  src="screenshots/image14.jpg"></A><BR>       
       <A href="screenshots/image14.jpg">screenshots/image6.jpg</A>
    </TD>  
</TR>
</TBODY>
</TABLE>





<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 
<H2>Demonstration of Implemented Features</H2>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 

<!------------------------------------------------------------------------> 
<H2><A name="texture"> 3D perspective viewing and objects </A></H2>
<!------------------------------------------------------------------------> 


Our game is built in 3D using R3Scene.cpp as a base, OpenGL libraries to set up the basic game properties, and raytrace.cpp and particle.cpp for various functions related to collisions and particle creations. The game is created primarily in a single XY plane, and the movement is entirely 2D, so the game genre can be considered 2.5D. This is very useful for the game editor (so a simple click can create an object only in the XY plane, rather than having to specify distance from camera). The camera that follows the player shifts based on the speed of the player.  

<hr><br>

<!------------------------------------------------------------------------> 
<H2><A name="texture"> Lighting and smooth shading.</A></H2>
<!------------------------------------------------------------------------> 

<!-- <TABLE>
<TBODY>
  <TR>
    <TD vAlign=top align=middle>
       <A href="output/particlemotion.gif"><IMG  width=256 src="output/particlemotion.gif"></A><BR>       
       <A href="input/particlemotion.scn">particlemotion.scn</A>
    </TD>
    <TD vAlign=top align=middle>
       <A href="output/effectofgravity.gif"><IMG width=256  src="output/effectofgravity.gif"></A><BR>       
       <A href="input/effectofgravity.scn">effectofgravity.scn</A>
    </TD>  
    <TD vAlign=top align=middle>
       <A href="output/effectofdrag.gif"><IMG width=256 src="output/effectofdrag.gif"></A><BR>       
       <A href="input/effectofdrag.scn">effectofdrag.scn</A>
    </TD>  
</TR>
</TBODY>
</TABLE> -->
We used OpenGL to implement lighting and smooth shading. Directional lights in the scn files light up the scene. The materials for all of the objects are defined in the scene. 
<hr><br>

<!------------------------------------------------------------------------> 
<H2><A name="texture"> User Input</A></H2>
<!------------------------------------------------------------------------> 

<!-- <TABLE>
<TBODY>
  <TR>
    <TD vAlign=top align=middle>
       <A href="output/particlemotion.gif"><IMG  width=256 src="output/particlemotion.gif"></A><BR>       
       <A href="input/particlemotion.scn">particlemotion.scn</A>
    </TD>
    <TD vAlign=top align=middle>
       <A href="output/effectofgravity.gif"><IMG width=256  src="output/effectofgravity.gif"></A><BR>       
       <A href="input/effectofgravity.scn">effectofgravity.scn</A>
    </TD>  
    <TD vAlign=top align=middle>
       <A href="output/effectofdrag.gif"><IMG width=256 src="output/effectofdrag.gif"></A><BR>       
       <A href="input/effectofdrag.scn">effectofdrag.scn</A>
    </TD>  
</TR>
</TBODY>
</TABLE> -->
The user controls the player character with three keys: "w", "a", and "d" mapped to up, left, and right.  Many other keys are used for auxiliary functions: "c" brings up the overview of the entire level, "m" toggles the minimap, "p" prints out a new scene file with the entire game in its current state, and the level editor allows the user to click anywhere in the scene and immediately create (or alter) objects. "r" restarts the level. 
<hr><br>

<!------------------------------------------------------------------------> 
<H2><A name="texture"> Computer Control over some element in the scene </A></H2>
<!------------------------------------------------------------------------> 

<!-- <TABLE>
<TBODY>
  <TR>
    <TD vAlign=top align=middle>
       <A href="output/particlemotion.gif"><IMG  width=256 src="output/particlemotion.gif"></A><BR>       
       <A href="input/particlemotion.scn">particlemotion.scn</A>
    </TD>
    <TD vAlign=top align=middle>
       <A href="output/effectofgravity.gif"><IMG width=256  src="output/effectofgravity.gif"></A><BR>       
       <A href="input/effectofgravity.scn">effectofgravity.scn</A>
    </TD>  
    <TD vAlign=top align=middle>
       <A href="output/effectofdrag.gif"><IMG width=256 src="output/effectofdrag.gif"></A><BR>       
       <A href="input/effectofdrag.scn">effectofdrag.scn</A>
    </TD>  
</TR>
</TBODY>
</TABLE> -->

We have an enemy AI for each of the enemies that is able to track the player's location and navigate towards him. Pairing this with a jumping ability and a lot of size, this creates significant challenges to the player navigating through the scene. We also create a fire particle source that creates particles that slowly fade from yellow to red. 

<hr><br>


<!------------------------------------------------------------------------> 
<H2><A name="texture"> Texture Mapping </A></H2>
<!------------------------------------------------------------------------> 

<!-- <TABLE>
<TBODY>
  <TR>
    <TD vAlign=top align=middle>
       <A href="output/particlemotion.gif"><IMG  width=256 src="output/particlemotion.gif"></A><BR>       
       <A href="input/particlemotion.scn">particlemotion.scn</A>
    </TD>
    <TD vAlign=top align=middle>
       <A href="output/effectofgravity.gif"><IMG width=256  src="output/effectofgravity.gif"></A><BR>       
       <A href="input/effectofgravity.scn">effectofgravity.scn</A>
    </TD>  
    <TD vAlign=top align=middle>
       <A href="output/effectofdrag.gif"><IMG width=256 src="output/effectofdrag.gif"></A><BR>       
       <A href="input/effectofdrag.scn">effectofdrag.scn</A>
    </TD>  
</TR>
</TBODY>
</TABLE> -->

We implemented texture mapping using OpenGL's libraries for texture maps. The texture file is passed
as an argument to the material in the scene file, which is then loaded onto the list of materials and applied onto an object through OpenGL. We texture mapped all players and enemies, coins, platforms, and the skybox behind the object (but only in one direction, to save memory). We also texture mapped the images in the sidebar of the level editor. 

<hr><br>

<!------------------------------------------------------------------------> 
<H2><A name="multipleviews"> Multiple Views </A></H2>
<!------------------------------------------------------------------------> 

<!-- <TABLE cellSpacing=0 cellPadding=10 border=0>
<TBODY>
  <TR>
    <TD vAlign=top align=middle>
       <A href="output/spheresource.gif"><IMG width=256 src="output/spheresource.gif"></A><BR>       
       <A href="input/spheresource.scn">spheresource.scn</A>
    </TD>
    <TD vAlign=top align=middle>

    </TD>
    <TD vAlign=top align=middle>
       <A href="output/meshsource.gif"><IMG width=256 src="output/meshsource.gif"></A><BR>       
       <A href="input/meshsource.scn">meshsource.scn</A>
    </TD>
  </TR>
</TBODY>
</TABLE>
 -->
We created a minimap on the bottom left corner that allows viewing of the entire map at the same time as the actual game, with a transparent background so it blends easily into the game. Furthermore, holding down "c" zooms the camera out to the entire map in fullscreen. 
<hr><br>

<!------------------------------------------------------------------------> 
<H2><A name="control"> On-Screen Control Panel </A></H2>
<!------------------------------------------------------------------------> 
<!-- 
<TABLE cellSpacing=0 cellPadding=10 border=0>
<TBODY>
  <TR>
    <TD vAlign=top align=middle>
       <A href="output/spheresink.gif"><IMG width=256 src="output/spheresink.gif"></A><BR>       
       <A href="input/spheresink.scn">spheresink.scn</A>
    </TD>
<!--     <TD vAlign=top align=middle>
       <A href="output/circlesink.gif"><IMG width=256 src="output/circlesink.gif"></A><BR>       
       <A href="input/circlesink.scn">circlesink.scn</A>
    </TD> -->
<!--     <TD vAlign=top align=middle>
       <A href="output/boxsink.gif"><IMG width=256 src="output/boxsink.gif"></A><BR>       
       <A href="input/boxsink.scn">boxsink.scn</A>
    </TD>
    <TD vAlign=top align=middle>
       <A href="output/linesink.gif"><IMG width=256 src="output/linesink.gif"></A><BR>       
       <A href="input/linesink.scn">linesink.scn</A>
    </TD>
  <TR>
 </TR> -->
<!--     <TD vAlign=top align=middle>
       <A href="output/cylindersink.gif"><IMG width=256 src="output/cylindersink.gif"></A><BR>       
       <A href="input/cylindersink.scn">cylindersink.scn</A>
    </TD>
    <TD vAlign=top align=middle>
       <A href="output/conesink.gif"><IMG width=256 src="output/conesink.gif"></A><BR>       
       <A href="input/conesink.scn">conesink.scn</A>
    </TD>
    <TD vAlign=top align=middle>
       <A href="output/meshsink.gif"><IMG width=256 src="output/meshsink.gif"></A><BR>       
       <A href="input/meshsink.scn">meshsink.scn</A>
    </TD> -->
<!--   </TR>
</TBODY>
</TABLE> --> 


There are multiple examples of on-screen control panels in our game. In the top left is a counter of the number of coins collected by the player. In the bottom left is a minimap HUD that allows the user to see the entire map at a time. In the level editor, there is a sidebar on the right side that allows you to choose between camera mode (for viewing), adding platforms, add enemies, 
<hr><br>

<!------------------------------------------------------------------------> 
<H2><A name="procedural"> Procedural and physically-based modeling </A></H2>
<!------------------------------------------------------------------------> 

We created particle sources that emanated particles that changed colors slowly from yellow to red, imitating fire as a decorative addition. We also created red GL_QUADs that fly out from the player in all directions once the player dies, simulating an explosion. 
<HR><BR>


<!------------------------------------------------------------------------> 
<H2><A name="particleinteraction">Collision Detection </A></H2>
<!------------------------------------------------------------------------> 

We implemented collision detection between the player and other objects, and between the enemies and other objects (including other enemies and players). All of the collisions are axis-aligned. We also detect the direction of collisions (e.g. the player hits the enemy from the top), so that the player can kill and be killed by enemies depending on the direction of approach. Collision detection handles impossible setups smoothly and without crashing (e.g. crushing blocks between platforms), though sometimes it may jump the block out in the wrong direction given an impossible situation. Collisions between moving objects are handled realistically: a jumping enemy on a jumping enemy will jump higher from the extra velocity. 

<HR><BR>

<!------------------------------------------------------------------------> 
<H2><A name="dynamics"> Simulated Dynamics </A></H2>
<!------------------------------------------------------------------------> 


We have several examples of simulated dynamics. The player moves left and right on a platform, accelerating to a max_speed, and jumps by adding a force in the upwards direction. The movement left and right is slowed by drag, so the player continues moving for a short time after holding left/right. Platforms accelerate away from their start points and slow towards the end points, and the player on a platform has equivalent acceleration (e.g. moving up quickly increases jump velocity). Enemies follow the same simulated dynamics and also have the ability to jump and follow the player wherever the player is, reacting dynamically to their environment. Fire particles go from yellow to red as they age.  
<HR><BR>


<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 
<H2><A name='sound'>Sound</A></H2>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 


We used the irrKlang library for sound, allowing wav files to be played whenever events occur. The user can select music to be played in the background in the scene file. Sounds are played during jumps, and the sound of enemy's jumps attenuates with the player's distance from the enemies. Collecting coins creates a Mario sound, death sounds are played when the player dies, and success sounds are played when a level is completed. 
<HR><BR>

<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 
<H2><A name='sound'>Level Editor</A></H2>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 

The level editor is a large proportion of our project. With the level editor, you can change the camera view to anywhere in the scene. You can select any block, platform, or enemy, and drag it to a different place using the mouse, as well as delete the object. You can select an enemy and make it bigger or smaller by dragging the mouse. You can also spawn enemies, stationary platforms, coins, or fire sources using only a mouse click at the location desired. All of these features are implemented using a convenient sidebar picker with 7 options, where you can simply click the option you want and then immediately begin creating objects. 
<BR>
Furthermore, we can immediately start the game at any time in the level editor by immediately pressing "n" to initiate the game. The current state of the game be printed at any time (requiring a lot of reverse engineering of Read()) by simply pressing "p", in which case a new .scn file is automatically generated with the location of all players, objects, fire sources, enemies, platforms, and coins, along with their associated textures. 
<HR><BR>
<BR>
<BR>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 
<H2><A name='sound'> Other Features</A></H2>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 

<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 
<H2><A name='sound'>SCN File</A></H2>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 

Our scn file has support for practically every variable in the game. Blocks and platforms can be created with any texture, and platforms move between two specified locations. The player and enemies can be created with any box shape at any location, with a given max speed. Enemies can be initialized to jump, to jump with any velocity, to follow the player, or to initially move in some direction (if it does not follow the player). Coins and fire sources can be created anywhere. The skybox image, minimum scene height, soundtrack to the scene, and next level after completion can all be set in the scene. 
<HR><BR>
<BR>


<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 
<H2><A name='sound'> Scene Printing </A></H2>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 

All relevant scene information can be immediately printed back out to another file ../levels/output.scn using "p", by reverse engineering Read() for all relevant objects. This is very useful for saving game state and creating new levels using the level editor. 

<HR><BR>
<BR>

<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 
<H2><A name='sound'> Multiple Level Integration </A></H2>
<!------------------------------------------------------------------------> 
<!------------------------------------------------------------------------> 

In our game, completing a level will immediately jump you to the next level as specified in the scn file. This transition is smooth, and can be used by level designers to build an integrated massive game (e.g. transitioning betwen many levels or dungeons). 

<HR><BR>
<BR>




<!------------------------------------------------------------------------> 

</BODY>
</HTML>
