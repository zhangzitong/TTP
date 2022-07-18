package utils;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class ConfigHelper {

    public static String getProperty(String name){
        //get datasets folder path
        Properties prop = new Properties();
        InputStream input;
        String ttpData = null;
        try {
            input = new FileInputStream("config.properties");
            // load a properties file
            prop.load(input);

            // get the property value and print it out
            ttpData = prop.getProperty(name);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ttpData;
    }

}
